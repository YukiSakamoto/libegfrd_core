
// Sakamoto EPDP Sample

#include <egfrd_core/config.h>
#include <stdexcept>
#include <vector>
#include <string>
#include <numeric>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <cstdlib>

// epdp headers
#include <utils/range.hpp>
#include <World.hpp>
#include <ParticleModel.hpp>
#include <SpeciesType.hpp>
#include <SpeciesTypeID.hpp>
#include <CuboidalRegion.hpp>
#include <NetworkRules.hpp>
#include <ReactionRule.hpp>
#include <EGFRDSimulator.hpp>
#include <GSLRandomNumberGenerator.hpp>

typedef double Real;

typedef ::World< ::CyclicWorldTraits<Real, Real> > world_type;

void throw_in_particles(void)
{

}

double distance_sq(world_type::position_type p1, world_type::position_type p2)
{
    double dsq = 0.0;
    world_type::position_type sq(gsl_pow_2(p2[0] - p1[0]), gsl_pow_2(p2[1] - p1[1]), gsl_pow_2(p2[2] - p2[2]));
    return std::accumulate(sq.begin(), sq.end(), 0.0);
}

class TemporaryParticleContainer {
public:
    typedef std::vector<std::pair< boost::shared_ptr< ::SpeciesType>, world_type::position_type> >  particle_position_container;
    TemporaryParticleContainer(void) {;}

    /*
    particle_position_container 
    list_particles_within_radius(boost::shared_ptr< ::SpeciesType> st, world_type::position_type &pos)
    {
        particle_position_container ret;
        for(particle_position_container::iterator it = container_.begin(); it != container_.end(); it++) {
            if (distance_sq(it->second, pos) < gsl_pow_2(std::atof((*st)["radius"].c_str() ) + std::atof((it->first)["radius"].c_str()))) 
            { 
                ret.push_back(*it); 
            }
        }
        return ret;
    }
    */

private:
    particle_position_container container_;
};


int main(int argc, char **argv)
{
    typedef ::World< ::CyclicWorldTraits<Real, Real> > world_type;
    typedef ::ParticleModel particle_model_type;
    typedef EGFRDSimulator< ::EGFRDSimulatorTraitsBase<world_type> > simulator_type;

    typedef ::CuboidalRegion<simulator_type::traits_type> cuboidal_region_type;
    typedef world_type::traits_type::structure_id_type structure_id_type;

    const Real world_size(1e-6);
    const Integer matrix_size(3);

    boost::shared_ptr<world_type> world(new world_type(world_size, matrix_size));
    world_type::position_type edge_length(world_size, world_size, world_size);
    world_type::position_type pos(world_size / 2, world_size / 2, world_size / 2);
    world->add_structure( boost::shared_ptr<cuboidal_region_type>(
                new cuboidal_region_type("world", cuboidal_region_type::shape_type(pos, pos))));
    boost::shared_ptr<GSLRandomNumberGenerator> rng(new GSLRandomNumberGenerator());
    particle_model_type model;
    rng->seed( (unsigned long int)0 );

    // add ::SpeciesType to ::ParticleModel
    boost::shared_ptr< ::SpeciesType> st(new ::SpeciesType());
    (*st)["name"] = std::string("A");
    (*st)["D"] = std::string("1e-12");
    (*st)["radius"] = std::string("2.5e-9");
    model.add_species_type(st);

    //add ::SpeciesInfo to ::World 
    const std::string &structure_id((*st)["structure"]);
    world->add_species( world_type::traits_type::species_type(
                st->id(), 
                boost::lexical_cast<world_type::traits_type::D_type>( (*st)["D"] ),
                boost::lexical_cast<world_type::length_type>( (*st)["radius"] ),
                boost::lexical_cast<structure_id_type>( structure_id.empty() ? "world" : structure_id )));

    int number_of_particles_A(10);
    for (int cnt = 0; cnt < number_of_particles_A; cnt++) {
        // add particles at random.
        world_type::position_type particle_pos( 
                rng->uniform(0.0, edge_length[0]), 
                rng->uniform(0.0, edge_length[1]), 
                rng->uniform(0.0, edge_length[2]) 
                );
        world->new_particle(st->id(), particle_pos);
    }


    return 0;
}
