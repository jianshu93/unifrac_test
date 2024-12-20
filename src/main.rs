use anyhow::Result;
use clap::{Arg, Command};
use phylotree::tree::Tree;
use std::path::Path;
use itertools::Itertools;
use unifrac::{
    compute::compute_unifrac_for_pair,
    io::{read_sample_table, write_matrix},
};

fn main() -> Result<()> {
    // Initialize logger
    println!("\n ************** initializing logger *****************\n");
    env_logger::Builder::from_default_env().init();
    let matches = Command::new("Unweighted_UniFrac")
        .version("0.1.0")
        .about("Fast Unweighted UniFrac")
        .arg(
            Arg::new("tree")
                .short('t')
                .long("tree")
                .value_name("TREE_FILE")
                .help("Input newick format tree file")
                .required(true),
        )
        .arg(
            Arg::new("table")
                .short('i')
                .long("input")
                .value_name("TABLE_FILE")
                .help("Input tab-delimited sample-feature table")
                .required(true),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("OUTPUT_FILE")
                .help("Output file for distance matrix")
                .required(true),
        )
        .get_matches();

    let tree_file = matches.get_one::<String>("tree").unwrap();
    let table_file = matches.get_one::<String>("table").unwrap();
    let output_file = matches.get_one::<String>("output").unwrap();

    // Read the tree
    let tree = Tree::from_file(Path::new(tree_file))?;

    // Read the sample-feature table
    let (taxa_order, sample_names, presence_matrix) = read_sample_table(table_file)?;
    assert!(
        presence_matrix.iter().map(|row| row.len()).all_equal(),
        "rows of the presence matrix are not all the same size..."
    );
    let n_samples = sample_names.len();

    // Compute distance matrix: n_samples x n_samples
    let mut dist_matrix = vec![0.0; n_samples * n_samples];

    for i in 0..n_samples {
        dist_matrix[i * n_samples + i] = 0.0; // distance to itself = 0
        for j in i + 1..n_samples {
            let uni = compute_unifrac_for_pair(&tree, &taxa_order, &presence_matrix, i, j)?;
            dist_matrix[i * n_samples + j] = uni;
            dist_matrix[j * n_samples + i] = uni; // symmetric
        }
    }

    // Write output matrix
    write_matrix(&sample_names, &dist_matrix, n_samples, output_file)?;

    Ok(())
}
