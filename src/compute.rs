use anyhow::{Context, Result};
use ndarray::{Array1, Array2, Zip};
use phylotree::tree::Tree;
use rayon::prelude::*;

/// Compute UniFrac for a given pair of samples i,j
pub fn compute_unifrac_for_pair(
    tree: &Tree,
    taxa_order: &[String],
    presence_matrix: &[Vec<f64>],
    i: usize,
    j: usize,
) -> Result<f64> {
    // Determine which taxa are present in either sample i or j
    let mut present_taxa = Vec::new();
    for (t_idx, taxon) in taxa_order.iter().enumerate() {
        let val_i = presence_matrix[t_idx][i];
        let val_j = presence_matrix[t_idx][j];
        if val_i > 0.0 || val_j > 0.0 {
            present_taxa.push(taxon.clone());
        }
    }

    let mut sub_tree = tree.clone();
    // prune taxa not in present_taxa
    {
        let leaves = sub_tree.get_leaves();
        let set: std::collections::HashSet<_> = present_taxa.iter().cloned().collect();
        for l in leaves {
            let name = sub_tree.get(&l).unwrap().name.clone().unwrap();
            if !set.contains(&name) {
                sub_tree.prune(&l).context("Prune failed")?;
                sub_tree.compress()?; // We also need to compress before pruning other leaves
            }
        }
    }

    let leaves = sub_tree.get_leaves();
    let mut leaf_order = vec![0; sub_tree.size()];
    let mut leaf_names = Vec::new();
    for (l_ord, l_idx) in leaves.into_iter().enumerate() {
        leaf_order[l_idx] = l_ord;
        leaf_names.push(sub_tree.get(&l_idx).unwrap().name.clone().unwrap());
    }

    let (mat_b, brlens) = construct_b(&sub_tree, &leaf_order)?;

    let p_a = get_sample_vec(&mat_b, presence_matrix, taxa_order, &leaf_names, i)?;
    let p_b = get_sample_vec(&mat_b, presence_matrix, taxa_order, &leaf_names, j)?;

    let sum_shared = parallel_elementwise_sum(&p_a, &p_b, &brlens);
    let l_total = brlens.sum();
    let unifrac = 1.0 - (sum_shared / l_total);

    Ok(unifrac)
}

/// Construct B and brlens
pub fn construct_b(tree: &Tree, leaf_order: &[usize]) -> Result<(Array2<u8>, Array1<f64>)> {
    let n_tips = tree.n_leaves();
    let n_branches = tree.size();
    let root = tree.get_root()?;

    let mut mat_b = Array2::<u8>::zeros((n_branches, n_tips));
    let mut brlens = Array1::zeros(n_branches);
    for idx in tree.postorder(&root)? {
        let node = tree.get(&idx)?;
        brlens[idx] = node.parent_edge.unwrap_or_default();
        if node.is_tip() {
            let t_ord = leaf_order[idx];
            mat_b[(idx, t_ord)] = 1;
        } else {
            for c in node.children.iter() {
                let merged = &mat_b.row(idx) + &mat_b.row(*c);
                mat_b.row_mut(idx).assign(&merged);
            }
        }
    }

    Ok((mat_b, brlens))
}

/// Construct p_a (or p_b) for a given sample index
pub fn get_sample_vec(
    mat: &Array2<u8>,
    presence_matrix: &[Vec<f64>],
    taxa_order: &[String],
    leaf_names: &[String],
    sample_idx: usize,
) -> Result<Array1<f64>> {
    let s = mat.shape();
    let mut p: Array1<f64> = Array1::zeros(s[0]);

    // For each leaf_name, find its taxon index in taxa_order, check presence in sample_idx
    for (col, lname) in leaf_names.iter().enumerate() {
        let t_idx = taxa_order.iter().position(|x| x == lname).unwrap();
        let val = presence_matrix[t_idx][sample_idx];
        if val > 0.0 {
            // Convert u8 to f64 before addition
            p = &p + &mat.column(col).mapv(|x| x as f64);
        }
    }

    // clamp to 0,1
    Ok(p.mapv(|v: f64| if v > 0.0 { 1.0 } else { 0.0 }))
}

/// Parallelize the element-wise multiply and sum (p_a * p_b * brlens)
pub fn parallel_elementwise_sum(p_a: &Array1<f64>, p_b: &Array1<f64>, brlens: &Array1<f64>) -> f64 {
    let data = p_a
        .iter()
        .zip(p_b.iter())
        .zip(brlens.iter())
        .map(|((pa, pb), b)| (*pa, *pb, *b))
        .collect::<Vec<_>>();

    data.par_iter().map(|(pa, pb, b)| pa * pb * b).sum()
}

/// Parallelize the element-wise multiply and sum (p_a * p_b * brlens)
pub fn parallel_elementwise_sum2(
    p_a: &Array1<f64>,
    p_b: &Array1<f64>,
    brlens: &Array1<f64>,
) -> f64 {
    Zip::from(p_a)
        .and(p_b)
        .and(brlens)
        .par_map_collect(|a, b, l| a * b * l) // Also collecting here...
        .sum()
}

pub fn vectorized_elementwise_sum(
    p_a: &Array1<f64>,
    p_b: &Array1<f64>,
    brlens: &Array1<f64>,
) -> f64 {
    (p_a * p_b * brlens).sum()
}
