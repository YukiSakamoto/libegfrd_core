
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

    void add( boost::shared_ptr< ::SpeciesType> st, world_type::position_type pos)
    {
        this->container_.push_back( particle_position_container::value_type(st, pos));
    }

    particle_position_container 
    list_particles_within_radius(boost::shared_ptr< ::SpeciesType> st, world_type::position_type &pos)
    {
        particle_position_container ret;
        for(particle_position_container::iterator it = container_.begin(); it != container_.end(); it++) {
            /*
            if (distance_sq(it->second, pos) < gsl_pow_2(std::atof((*st)["radius"].c_str() ) + std::atof((it->first)["radius"].c_str()))) 
            { 
                ret.push_back(*it); 
            }
            */
            double radius_new( atof(((*st)["radius"]).c_str()) );
            double radius_st(  atof(((*(it->first))["radius"]).c_str()) );
            if (distance_sq(it->second, pos) < gsl_pow_2(radius_new) ) {
                ret.push_back( *it );
            }
        }
        return ret;
    }

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
    const Real volume( world_size * world_size * world_size);

    const Integer N(100);
    const Real kd(0.1), U(0.5);
    const Real ka(kd * volume * (1 - U) / (U * U * N));
    const Real k2(ka), k1(kd);

    boost::shared_ptr<world_type> world(new world_type(world_size, matrix_size));
    world_type::position_type edge_length(world_size, world_size, world_size);
    world_type::position_type pos(world_size / 2, world_size / 2, world_size / 2);
    world->add_structure( boost::shared_ptr<cuboidal_region_type>(
                new cuboidal_region_type("world", cuboidal_region_type::shape_type(pos, pos))));
    boost::shared_ptr<GSLRandomNumberGenerator> rng(new GSLRandomNumberGenerator());
    particle_model_type model;
    rng->seed( (unsigned long int)0 );

    // add ::SpeciesType to ::ParticleModel
    boost::shared_ptr< ::SpeciesType> st1(new ::SpeciesType());
    {
        (*st1)["name"] = std::string("A");
        (*st1)["D"] = std::string("1e-12");
        (*st1)["radius"] = std::string("2.5e-9");
    }
    model.add_species_type(st1);

    boost::shared_ptr< ::SpeciesType> st2(new ::SpeciesType());
    {
        (*st2)["name"] = std::string("A");
        (*st2)["D"] = std::string("1e-12");
        (*st2)["radius"] = std::string("2.5e-9");
    }
    model.add_species_type(st2);

    boost::shared_ptr< ::SpeciesType> st3(new ::SpeciesType());
    {
        (*st3)["name"] = std::string("A");
        (*st3)["D"] = std::string("1e-12");
        (*st3)["radius"] = std::string("2.5e-9");
    }
    model.add_species_type(st3);

    // A -> B + C   k1
    std::vector< ::SpeciesTypeID> products;
    products.push_back(st2->id());
    products.push_back(st3->id());
    model.network_rules().add_reaction_rule( new_reaction_rule(st1->id(), products, k1) );

    // B + C -> A   k2
    products.clear();
    products.push_back(st1->id());
    model.network_rules().add_reaction_rule( new_reaction_rule(st2->id(), st3->id(), products, k2) );

    {   // world::set_all_repusive() equality section
        BOOST_FOREACH( boost::shared_ptr< ::SpeciesType> temp_st1, model.get_species_types()) {
            BOOST_FOREACH( boost::shared_ptr< ::SpeciesType> temp_st2, model.get_species_types()) {
                boost::scoped_ptr< ::NetworkRules::reaction_rule_generator> gen( model.network_rules().query_reaction_rule( temp_st1->id(), temp_st2->id()));
                if (!gen) {
                    const::std::vector< ::SpeciesTypeID> products;
                    model.network_rules().add_reaction_rule( ::new_reaction_rule(temp_st1->id()(), temp_st2->id(), products, 0.0) );
                }
            }
        }
    }

    //add ::SpeciesInfo to ::World 
    const std::string &structure_id((*st1)["structure"]);
    world->add_species( world_type::traits_type::species_type(
                st1->id(), 
                boost::lexical_cast<world_type::traits_type::D_type>( (*st1)["D"] ),
                boost::lexical_cast<world_type::length_type>( (*st1)["radius"] ),
                boost::lexical_cast<structure_id_type>( structure_id.empty() ? "world" : structure_id )));

    int number_of_particles_A(N);
    TemporaryParticleContainer container;
    for (int cnt = 0; cnt < number_of_particles_A; cnt++) {
        // add particles at random.
        for(;;) {
            world_type::position_type particle_pos( rng->uniform(0.0, edge_length[0]), rng->uniform(0.0, edge_length[1]), rng->uniform(0.0, edge_length[2]) );
            if (container.list_particles_within_radius(st1, particle_pos).size() == 0) {
                //std::cout << "(" << particle_pos[0] << particle_pos[1] << particle_pos[2] << ")" << std::endl;
                container.add(st1, particle_pos);
                world->new_particle(st1->id(), particle_pos);
                break;
            }
        }
    }

    return 0;
}
