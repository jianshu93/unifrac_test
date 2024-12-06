use anyhow::{Context, Result};
use ndarray::{Array1, Array2};
use phylotree::tree::Tree;
use std::{
    env::args,
    path::Path,
    time::{Duration, Instant},
};

fn main() -> Result<()> {
    let fname = args()
        .nth(1)
        .context("Must give a newick tree file as input")?;

    let tree = Tree::from_file(Path::new(&fname))?;

    // This could be an iterator method on the tree.nodes vector
    // if it was `pub`
    let mut leaf_order = vec![0; tree.size()];
    let mut labs = vec![];
    for (l_ord, l_idx) in tree.get_leaves().into_iter().enumerate() {
        leaf_order[l_idx] = l_ord;
        labs.push(tree.get(&l_idx).unwrap().name.clone().unwrap());
    }

    // Cosntruct the matrix B
    let (mat_b, brlens) = construct_b(&tree, &leaf_order[..])?;

    eprintln!("Tree: {}", tree.to_newick()?);
    eprintln!("Leaf order = {:?}", labs); // Order of the leaves/B's columns
    eprintln!("B =\n{}", mat_b);
    eprintln!("b = {}", brlens);

    let s_a = [0, 1, 2, 3]; // T1,T2,T3,T4
    let s_b = [1, 2]; // T2, T3

    let p_a = get_sample_vec(&mat_b, &s_a[..])?;
    let p_b = get_sample_vec(&mat_b, &s_b[..])?;

    eprintln!("\nP_A:{p_a}");
    eprintln!("P_B:{p_b}");

    // Computation with the matrix product
    let br = Array2::eye(brlens.shape()[0]) * &brlens; // Diagonal matrix with branch lengths

    let mut times = vec![];
    for _ in 0..1000 {
        let start = Instant::now();
        let _ = p_a.t().dot(&br).dot(&p_b);
        times.push(start.elapsed());
    }
    let d = times.iter().sum::<Duration>() / 1000;

    println!("\nComputing unifrac:");

    let f = p_a.t().dot(&br).dot(&p_b);

    let mut times = vec![];
    for _ in 0..1000 {
        let p_a_2 = p_a.clone();
        let start = Instant::now();
        let _ = (p_a_2 * &p_b * &brlens).sum();
        times.push(start.elapsed());
    }
    let d_2 = times.iter().sum::<Duration>() / 1000;

    // Computation without matrix product
    let f2 = (p_a * &p_b * &brlens).sum();
    eprintln!("P_A.T @ Br @ P_B   = {f}\t({d:?})");
    eprintln!("Sum(P_A * P_B * b) = {f2}\t({d_2:?})");

    let unifrac = 1. - f2 / brlens.sum();
    eprintln!("D_unifrac(S_A,S_B) = {unifrac:.3}");

    Ok(())
}

/// Construct the b matrix of branches that are between the root and a given taxa,
/// As well as the corresponding vector of branch lengths.
/// This returns:
///  - B an (n_branches X n_tips) binary matrix where a 1 in cell (i,j) indicates that
///    branch i is in the path between the root and leaf j.
///  - b a (n_branches X 1) vector containing the branch lenghts in the same order as B
fn construct_b(tree: &Tree, leaf_order: &[usize]) -> Result<(Array2<u8>, Array1<f64>)> {
    let n_tips = tree.n_leaves();
    let n_branches = tree.size();
    let root = tree.get_root()?;

    let mut mat_b = Array2::<u8>::zeros((n_branches, n_tips));
    let mut brlens = Array1::zeros(n_branches);
    for idx in tree.postorder(&root)? {
        let node = tree.get(&idx)?;
        brlens[idx] = node.parent_edge.unwrap_or_default();
        if node.is_tip() {
            // Set the tip bit to 1
            let t_ord = leaf_order[idx];
            mat_b[(idx, t_ord)] = 1;
        } else {
            // Add the descendants
            for c in node.children.iter() {
                let merged = &mat_b.row(idx) + &mat_b.row(*c);
                mat_b.row_mut(idx).assign(&merged);
            }
        }
    }

    Ok((mat_b, brlens))
}

/// Constructs the binary vector of branches contained in a sample
/// (i.e. the union of branches between the root and each tip in the sample)
fn get_sample_vec(mat: &Array2<u8>, sample: &[usize]) -> Result<Array1<f64>> {
    let s = mat.shape();
    let mut p_a = Array1::zeros(s[0]);
    for i in sample {
        p_a += &mat.column(*i);
    }

    Ok(p_a.clamp(0, 1).mapv(|v| v as f64))
}
