use clap::Parser;
use forrustts::genetics::GeneticMap;
use forrustts::genetics::GeneticMapBuilder;
use forrustts::genetics::PoissonCrossover;

use design_playground::genome_array::*;
use design_playground::SimParams;

fn main() {
    let params = SimParams::parse();
    let genome_length = params.genome_length.unwrap();
    let poisson =
        vec![
            PoissonCrossover::new(0, genome_length, params.recrate * (genome_length as f64))
                .unwrap(),
        ];
    let builder = GeneticMapBuilder::default().extend_poisson(&poisson);
    let genetic_map = GeneticMap::new_from_builder(builder).unwrap();
    let pop = evolve_pop(params, genetic_map).unwrap();
    println!(
        "{} {}",
        pop.num_segregating_mutations(),
        pop.sum_extant_genome_sizes()
    );
}
