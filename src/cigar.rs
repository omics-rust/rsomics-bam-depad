use crate::bam_layout::{
    CIGAR_D, CIGAR_H, CIGAR_I, CIGAR_M, CIGAR_P, CIGAR_S, FIXED_HEAD, L_READ_NAME, N_CIGAR,
    cigar_gen,
};

/// Decode the original CIGAR ops from a BAM payload.
pub(crate) fn decode_cigar(record: &[u8]) -> Vec<(u8, u32)> {
    let name_len = usize::from(record[L_READ_NAME]);
    let n_cigar = usize::from(u16::from_le_bytes([record[N_CIGAR], record[N_CIGAR + 1]]));
    let cigar_start = FIXED_HEAD + name_len;
    (0..n_cigar)
        .map(|ci| {
            let raw = u32::from_le_bytes(
                record[cigar_start + ci * 4..cigar_start + ci * 4 + 4]
                    .try_into()
                    .unwrap(),
            );
            ((raw & 0xf) as u8, raw >> 4)
        })
        .collect()
}

/// H and S clip lengths at the 5' end.
pub(crate) fn leading_clips(cigar: &[(u8, u32)]) -> (u32, u32) {
    let mut h = 0u32;
    let mut s = 0u32;
    let mut i = 0usize;
    if !cigar.is_empty() && cigar[0].0 == CIGAR_H {
        h = cigar[0].1;
        i = 1;
    }
    if i < cigar.len() && cigar[i].0 == CIGAR_S {
        s = cigar[i].1;
    }
    (h, s)
}

/// S and H clip lengths at the 3' end.
pub(crate) fn trailing_clips(cigar: &[(u8, u32)]) -> (u32, u32) {
    let n = cigar.len();
    if n == 0 {
        return (0, 0);
    }
    let mut h = 0u32;
    let mut s = 0u32;
    let mut i = n;
    if cigar[n - 1].0 == CIGAR_H {
        h = cigar[n - 1].1;
        i -= 1;
    }
    if i > 0 && cigar[i - 1].0 == CIGAR_S {
        s = cigar[i - 1].1;
    }
    (s, h)
}

/// Rebuild the CIGAR for a query read after depadding.
///
/// Mirrors `bam_pad2unpad`'s CIGAR-reconstruction loop (padding.c).
pub(crate) fn rebuild_cigar(record: &[u8], q: &[u8], ref_seq: &[u8], pos: usize) -> Vec<u32> {
    let orig_cigar = decode_cigar(record);
    let mut out: Vec<u32> = Vec::with_capacity(q.len() + 4);

    let (lead_h, lead_s) = leading_clips(&orig_cigar);
    if lead_h > 0 {
        out.push(cigar_gen(lead_h, CIGAR_H));
        if lead_s > 0 {
            out.push(cigar_gen(lead_s, CIGAR_S));
        }
    } else if lead_s > 0 {
        out.push(cigar_gen(lead_s, CIGAR_S));
    }

    if q.is_empty() {
        let (trail_s, trail_h) = trailing_clips(&orig_cigar);
        if trail_s > 0 {
            out.push(cigar_gen(trail_s, CIGAR_S));
        }
        if trail_h > 0 {
            out.push(cigar_gen(trail_h, CIGAR_H));
        }
        return out;
    }

    let mut ops: Vec<u8> = Vec::with_capacity(q.len());
    for (i, &qb) in q.iter().enumerate() {
        let k = pos + i;
        let rb = ref_seq[k];
        ops.push(match (qb != 0, rb != 0) {
            (true, true) => CIGAR_M,
            (true, false) => CIGAR_I,
            (false, true) => CIGAR_D,
            (false, false) => CIGAR_P,
        });
    }

    // Count leading pad columns immediately before `pos`. Mirrors the
    // `include any pads if starts with an insert` block in padding.c.
    let lead_pad_k: u32 = if matches!(ops[0], CIGAR_I | CIGAR_P) {
        let mut scan = pos;
        let mut k = 0u32;
        while scan > 0 && ref_seq[scan - 1] == 0 {
            k += 1;
            scan -= 1;
        }
        k
    } else {
        0
    };

    if lead_pad_k > 0 && ops[0] == CIGAR_I {
        out.push(cigar_gen(lead_pad_k, CIGAR_P));
    }

    let start_len: u32 = if ops[0] == CIGAR_P { lead_pad_k + 1 } else { 1 };
    let mut cur_op = ops[0];
    let mut cur_len = start_len;

    for &op in &ops[1..] {
        if op == cur_op {
            cur_len += 1;
        } else {
            out.push(cigar_gen(cur_len, cur_op));
            cur_op = op;
            cur_len = 1;
        }
    }
    out.push(cigar_gen(cur_len, cur_op));

    remove_redundant_p(&mut out);

    let (trail_s, trail_h) = trailing_clips(&orig_cigar);
    if trail_h > 0 {
        if trail_s > 0 {
            out.push(cigar_gen(trail_s, CIGAR_S));
        }
        out.push(cigar_gen(trail_h, CIGAR_H));
    } else if trail_s > 0 {
        out.push(cigar_gen(trail_s, CIGAR_S));
    }

    out
}

/// Remove redundant P operators between M/D operators (padding.c post-processing pass).
///
/// `5M 2P 10M` → `15M`. If flanking ops differ only P is dropped without merging.
pub(crate) fn remove_redundant_p(cigar: &mut Vec<u32>) {
    let n = cigar.len();
    if n < 3 {
        return;
    }
    let mut i = 2;
    while i < cigar.len() {
        let pre_op = (cigar[i - 2] & 0xf) as u8;
        let mid_op = (cigar[i - 1] & 0xf) as u8;
        let post_op = (cigar[i] & 0xf) as u8;
        if mid_op == CIGAR_P
            && (pre_op == CIGAR_M || pre_op == CIGAR_D)
            && (post_op == CIGAR_M || post_op == CIGAR_D)
        {
            cigar[i - 1] = 0; // mark P as 0 (= 0M, filtered below)
            if pre_op == post_op {
                let combined = (cigar[i - 2] >> 4) + (cigar[i] >> 4);
                cigar[i] = cigar_gen(combined, post_op);
                cigar[i - 2] = 0;
            }
        }
        i += 1;
    }
    cigar.retain(|&op| op != 0);
}

/// Replace the CIGAR section of a raw BAM payload, returning a new byte vector.
///
/// Name/seq/qual/aux are copied verbatim; only CIGAR and `n_cigar` change.
pub(crate) fn replace_cigar(record: &[u8], new_cigar: &[u32]) -> Vec<u8> {
    let name_len = usize::from(record[L_READ_NAME]);
    let old_n_cigar = usize::from(u16::from_le_bytes([record[N_CIGAR], record[N_CIGAR + 1]]));
    let cigar_start = FIXED_HEAD + name_len;
    let after_old = cigar_start + old_n_cigar * 4;

    let mut out = Vec::with_capacity(
        FIXED_HEAD + name_len + new_cigar.len() * 4 + (record.len() - after_old),
    );
    out.extend_from_slice(&record[..FIXED_HEAD + name_len]);
    for &op in new_cigar {
        out.extend_from_slice(&op.to_le_bytes());
    }
    out.extend_from_slice(&record[after_old..]);
    let nc = new_cigar.len() as u16;
    out[N_CIGAR..N_CIGAR + 2].copy_from_slice(&nc.to_le_bytes());
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam_layout::cigar_gen;

    #[test]
    fn remove_redundant_p_merges() {
        // 5M 2P 10M → 15M
        let mut cigar = vec![
            cigar_gen(5, CIGAR_M),
            cigar_gen(2, CIGAR_P),
            cigar_gen(10, CIGAR_M),
        ];
        remove_redundant_p(&mut cigar);
        assert_eq!(cigar, vec![cigar_gen(15, CIGAR_M)]);
    }

    #[test]
    fn remove_redundant_p_md_removes_p() {
        // 5M 2P 3D: P between M and D is removed but ops differ so no merge.
        let mut cigar = vec![
            cigar_gen(5, CIGAR_M),
            cigar_gen(2, CIGAR_P),
            cigar_gen(3, CIGAR_D),
        ];
        remove_redundant_p(&mut cigar);
        assert_eq!(cigar, vec![cigar_gen(5, CIGAR_M), cigar_gen(3, CIGAR_D)]);
    }
}
