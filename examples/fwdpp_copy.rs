use clap::Parser;
use forrustts::genetics::GeneticMap;
use forrustts::genetics::GeneticMapBuilder;
use forrustts::genetics::PoissonCrossover;

use design_playground::fwdpp_copy::*;
use design_playground::SimParams;

fn main() {
    let params = SimParams::parse();
    let poisson = vec![PoissonCrossover::new(0, 1000000, params.recrate).unwrap()];
    let builder = GeneticMapBuilder::default().extend_poisson(&poisson);
    let genetic_map = GeneticMap::new_from_builder(builder).unwrap();
    let pop = evolve_pop_with_haplotypes(params, genetic_map).unwrap();
    println!(
        "{} {} {} {}",
        pop.num_segregating_mutations(),
        pop.sum_extant_genome_sizes(),
        pop.num_extant_genomes(),
        pop.total_num_genomes()
    );
}
