use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use csv::Writer;
use std::ffi::OsStr;
use std::fs;
use std::iter::zip;
use std::path::Path;

use orca::{molecule::Molecule, loader, assembly::{
    index_search, Bound, log_bound, addition_bound,
}};

pub fn dataset_bench(c: &mut Criterion) {
    // Define a new criterion benchmark group of dataset benchmarks.
    let mut group = c.benchmark_group("datasets");

    // Loop over all datasets of interest.
    for dataset in ["gdb13_1201", "gdb17_800"].iter() {
        // Load all molecules from the given dataset.
        let paths = fs::read_dir(Path::new("data").join(dataset)).unwrap();
        let mut mol_list: Vec<Molecule> = Vec::new();
        for path in paths {
            let name = path.unwrap().path();
            if name.extension().and_then(OsStr::to_str) != Some("mol") {
                continue;
            }
            mol_list.push(
                loader::parse_molfile_str(
                    &fs::read_to_string(name.clone())
                    .expect(&format!("Could not read file {name:?}"))
                ).expect(&format!("Failed to parse {name:?}")));
        }

        // For each of the bounds options, run the benchmark over all molecules
        // in this dataset.
        let bounds = [vec![], vec![Bound::Log(log_bound)], vec![Bound::Addition(addition_bound)]];
        let bound_strs = ["naive", "logbound", "addbound"];
        for (bound, bound_str) in zip(&bounds, &bound_strs) {
            group.bench_with_input(
                BenchmarkId::new(format!("{bound_str}"), &dataset),
                bound, |b, bound| {
                    b.iter(|| {
                        for mol in &mol_list {
                            index_search(&mol, &bound);
                        }
                    });
            });
        }
    }

    group.finish();
}

pub fn jossplot_bench(c: &mut Criterion) {
    // Define a new criterion benchmark group of dataset benchmarks.
    let mut group = c.benchmark_group("jossplot");

    // Set up CSV file for recording the number of duplicate isomorphic
    // subgraphs per molecule.
    let crit_path = Path::new("target").join("criterion").join("jossplot");
    fs::create_dir_all(&crit_path).unwrap();
    let mut csv = Writer::from_path(&crit_path.join("jossplot.csv")).unwrap();

    // Loop over all datasets of interest.
    for dataset in ["gdb13_1201", "gdb17_800"].iter() {
        // Iterate over each molecule file from the current dataset.
        let mut paths: Vec<_> =
            fs::read_dir(Path::new("data").join(dataset))
            .unwrap()
            .filter_map(|r| r.ok())
            .collect();
        paths.sort_by_key(|p| p.path());
        for path in paths {
            // If the file is a .mol, parse it as a Molecule.
            let name = path.path();
            if name.extension().and_then(OsStr::to_str) != Some("mol") {
                continue;
            }
            let mol = loader::parse_molfile_str(
                &fs::read_to_string(name.clone())
                .expect(&format!("Could not read file {name:?}"))
            ).expect(&format!("Failed to parse {name:?}"));

            // Make a string ID using this molecule's dataset and filename.
            let mol_fname = name.file_name().unwrap().to_str().unwrap();
            let id = format!("{dataset}-{mol_fname}");

            // Benchmark assembly index calculation for this molecule using the
            // different bound options.
            let bounds = [vec![], vec![Bound::Log(log_bound)], vec![Bound::Addition(addition_bound)]];
            let bound_strs = ["naive", "logbound", "addbound"];
            for (bound, bound_str) in zip(&bounds, &bound_strs) {
                group.bench_with_input(
                    BenchmarkId::new(format!("{bound_str}"), &id),
                    bound, |b, bound| {
                        b.iter(|| index_search(&mol, &bound));
                });
            }

            // Record this molecule's number of duplicate isomorphic subgraphs
            // to the CSV file.
            let dup_iso_subs = mol.matches().count().to_string();
            csv.write_record(&[id, dup_iso_subs]).unwrap();
        }
    }

    group.finish();
    csv.flush().unwrap();
}

criterion_group!(benches, dataset_bench, jossplot_bench);
criterion_main!(benches);
