use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use csv::Writer;
use std::ffi::OsStr;
use std::fs;
use std::iter::zip;

use orca::{molecule::Molecule, loader, assembly::{
    index_search, Bound, log_bound, addition_bound,
}};

pub fn dataset_bench(c: &mut Criterion) {
    // Define a new criterion benchmark group of dataset benchmarks.
    let mut group = c.benchmark_group("datasets");

    // Loop over all datasets of interest.
    for dataset in ["gdb13_1201", "gdb17_800"].iter() {
        // Load all molecules from the given dataset.
        let paths = fs::read_dir(&format!("data/{dataset}")).unwrap();
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
        let bound_strs = ["nobound", "logbound", "addbound"];
        for (bound, bound_str) in zip(&bounds, &bound_strs) {
            let id = format!("{dataset}-{bound_str}");
            group.bench_with_input(
                BenchmarkId::new("index_search", &id), bound, |b, bound| {
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
    let benchmark = "gdb17_800";
    fs::create_dir_all("./target/criterion").unwrap();
    let mut csv_wtr = Writer::from_path("./target/criterion/molecule_space.csv").unwrap();

    let paths = fs::read_dir(format!("./data/{benchmark}")).unwrap();
    let mut group = c.benchmark_group("plot");
    for path in paths {
        let file_path = path.unwrap();
        if file_path.file_type().unwrap().is_file() {
            let file_name = file_path.file_name();
            let file_name_clone = file_name.clone();
            let file = file_name_clone.into_string().unwrap();
            
            match file.split(".").nth(1).unwrap() {
                "mol" | "sdf" => {
                    let mol_name = file_name.to_str().unwrap().split(".").nth(0).unwrap();
                    let molfile = fs::read_to_string(&format!("./data/{benchmark}/{file}")).expect("Cannot read input file.");
                    let molecule = loader::parse_molfile_str(&molfile).expect("Cannot parse molfile.");
                    
                    // Run the three bounds for generating Assembly Index
                    group.bench_function(format!("{mol_name}_Naive"), |b| b.iter(|| index_search(&molecule, &[])));
                    group.bench_function(format!("{mol_name}_Log"),|b| b.iter(|| index_search(&molecule, &[Bound::Log(log_bound)])));
                    group.bench_function(format!("{mol_name}_Addition"),|b| b.iter(|| index_search(&molecule, &[Bound::Addition(addition_bound)])));
                    
                    // Run the search-space algo
                    let dup_subgraphs = molecule.matches().count().to_string();
                    csv_wtr.write_record(&[mol_name, &dup_subgraphs]).unwrap();
                }
                _ => {}
            }
        }
    }
    group.finish();
    csv_wtr.flush().unwrap();
}

criterion_group!(benches, dataset_bench, jossplot_bench);
criterion_main!(benches);
