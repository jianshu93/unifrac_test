use anyhow::{Context, Result};
use std::{
    fs::File,
    io::{BufRead, BufReader, Write},
};

/// Read the sample-feature table.
/// First line: ignore the first element, subsequent elements are sample names
/// Example:
/// Anything  SampleA  SampleB  SampleC
/// T1        10       0        5
/// T2        0        25       0
/// ...
///
/// Any value > 0 is converted to 1.0, else 0.0.
pub fn read_sample_table(filename: &str) -> Result<(Vec<String>, Vec<String>, Vec<Vec<f64>>)> {
    let f = File::open(filename)?;
    let mut lines = BufReader::new(f).lines();

    // First line: parse sample names
    let header = lines.next().context("No header in table")??;
    let mut hdr_split = header.split('\t');
    hdr_split.next(); // ignore the first element in the header line
    let sample_names: Vec<String> = hdr_split.map(|s| s.to_string()).collect();

    let mut taxa_order = Vec::new();
    let mut presence_matrix = Vec::new();

    for line in lines {
        let line = line?;
        let mut parts = line.split('\t');
        let taxon = parts.next().context("Taxon missing in a line")?.to_string();
        taxa_order.push(taxon);
        let values: Vec<f64> = parts
            .map(|x| {
                let val: f64 = x.parse().unwrap_or(0.0);
                if val > 0.0 {
                    1.0
                } else {
                    0.0
                } // ensure binary
            })
            .collect();
        presence_matrix.push(values);
    }

    Ok((taxa_order, sample_names, presence_matrix))
}

/// Write the resulting matrix to a file
pub fn write_matrix(
    sample_names: &[String],
    dist_matrix: &[f64],
    n: usize,
    output_file: &str,
) -> Result<()> {
    let mut file = File::create(output_file)?;
    // Write header
    // The first column header is sample names as well
    write!(file, "Sample")?;
    for sn in sample_names {
        write!(file, "\t{}", sn)?;
    }
    writeln!(file)?;

    for i in 0..n {
        write!(file, "{}", sample_names[i])?;
        for j in 0..n {
            write!(file, "\t{:.6}", dist_matrix[i * n + j])?;
        }
        writeln!(file)?;
    }

    Ok(())
}
