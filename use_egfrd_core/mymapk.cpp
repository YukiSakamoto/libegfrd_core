
// Sakamoto EPDP Sample

#include <stdexcept>
#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

// epdp headers
#include <utils/range.hpp>
#include <World.hpp>
#include <ParticleModel.hpp>
#include <SpeciesType.hpp>
#include <SpeciesTypeID.hpp>
#include <CuboidalRegion.hpp>
#include <NetworkRules.hpp>
#include <ReactionRule.hpp>
#include <CuboidalRegion.hpp>
#include <EGFRDSimulator.hpp>
#include <GSLRandomNumberGenerator.hpp>

typedef double Real;


int main(int argc, char **argv)
{
    typedef ::World< ::CyclicWorldTraits<Real, Real> > world_type;
    typedef ::ParticleModel particle_model_type;
    typedef EGFRDSimulator< ::EGFRDSimulatorTraitsBase<world_type> > simulator_type;

    typedef ::CuboidalRegion<simulator_type::traits_type> cuboidal_region_type;

    const Real world_size(1e-6);
    const Integer matrix_size(3);

    boost::shared_ptr<world_type> world(new world_type(world_size, matrix_size));
    world_type::position_type pos(world_size / 2, world_size / 2, world_size / 2);
    world->add_structure( boost::shared_ptr<cuboidal_region_type>(
                new cuboidal_region_type("world", cuboidal_region_type::shape_type(pos, pos))));
    boost::shared_ptr<GSLRandomNumberGenerator> rng(new GSLRandomNumberGenerator());
    rng->seed( (unsigned long int)0 );

    boost::shared_ptr< ::SpeciesType> st(new ::SpeciesType());
    (*st)["name"] = std::string("A");
    (*st)["D"] = std::string("1e-12");
    (*st)["radius"] = std::string("2.5e-9");

    return 0;
}
