use std::fs;

use orca::{assembly::{
    addition_bound, index, index_and_states, log_bound, search_space, Bound,
}, molecule::Molecule};
use csv::Writer;

use criterion::{criterion_group, criterion_main, Criterion};
use orca::loader;


pub fn plot_benchmark(c: &mut Criterion) {
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
                    group.bench_function(format!("{mol_name}_Naive"), |b| b.iter(|| index_and_states(&molecule, &[])));
                    group.bench_function(format!("{mol_name}_Log"),|b| b.iter(|| index_and_states(&molecule, &[Bound::Log(log_bound)])));
                    group.bench_function(format!("{mol_name}_Addition"),|b| b.iter(|| index_and_states(&molecule, &[Bound::Addition(addition_bound)])));
                    
                    // Run the search-space algo
                    csv_wtr.write_record(&[mol_name, &search_space(&molecule).to_string()]).unwrap();
                }
                _ => {}
            }
        }
    }
    group.finish();
    csv_wtr.flush().unwrap();
}

pub fn gdb13_benchmark(c: &mut Criterion) {
    let paths = fs::read_dir("data/gdb13_1201").unwrap(); //gdb13

    let mut molecules_list: Vec<Molecule> = Vec::new();
    for path in paths {
        let name = path.unwrap().path();
        let molfile = fs::read_to_string(name.clone()).expect("Cannot read file");
        let molecule = loader::parse_molfile_str(&molfile).expect("Cannot parse molecule");
        molecules_list.push(molecule);
    }
    c.bench_function("gdb13", |b| b.iter(|| 
        {
            for molecule in &molecules_list {
                index(&molecule);
            }
        }
    ));
}

pub fn gdb17_benchmark(c: &mut Criterion) {
    let paths = fs::read_dir("data/gdb17_800").unwrap(); //gdb13

    let mut molecules_list: Vec<Molecule> = Vec::new();
    for path in paths {
        let name = path.unwrap().path();
        let molfile = fs::read_to_string(name.clone()).expect("Cannot read file");
        let molecule = loader::parse_molfile_str(&molfile).expect("Cannot parse molecule");
        molecules_list.push(molecule);
    }
    c.bench_function("gdb17", |b| b.iter(|| 
        {
            for molecule in &molecules_list {
                index(&molecule);
            }
        }
    ));
}

criterion_group!(benches, gdb13_benchmark, gdb17_benchmark, plot_benchmark);
criterion_main!(benches);
