#ifndef NCPlugin_CrystallineExtinction_hh
#define NCPlugin_CrystallineExtinction_hh

#include "NCrystal/NCPluginBoilerplate.hh"//Common stuff (includes NCrystal
                                          //public API headers, sets up
                                          //namespaces and aliases)
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/extd_utils/NCPlaneProvider.hh"

namespace NCPluginNamespace {

  enum class Recipe { cls, std, lux };

  //We implement the actual physics model in this completely custom C++ helper
  //class. That decouples it from NCrystal interfaces (which is nice in case the
  //NCrystal API changes at some point), and it makes it easy to directly
  //instantiate and test the modelling implementation from standalone C++ code.
  //
  //We mark the class as MoveOnly, to make sure it doesn't get copied around by
  //accident (since it could easily end up having large data members).

  class CrystallineExtinction final : public NC::MoveOnly {
  public:

    //A few static helper functions which can extract relevant data from NCInfo
    //objects (the createFromInfo function will raise BadInput exceptions in
    //case of syntax errors in the @CUSTOM_ section data):

    static bool isApplicable( const NC::Info& );
    static CrystallineExtinction createFromInfo( const NC::Info&, NC::PlaneProvider * = nullptr );//will raise BadInput in case of syntax errors

    //To account for the block size and the mosaic spread, the theories of extinction
    //introduce primary and secondary corrections to the coherent elastic neutron
    //scattering cross section for randomly oriented powders. Dependent on whether
    //blocks separated by small-angle tilts but having the same relative orientations
    //are considered correlated, two models have been developed, the correlated block
    //model and the uncorrelated block model, respectively.
    //Ref: International Tables for Crystallography (2006). Vol. C, Chapter 6.4, pp. 609–616.
    //As a first investigation, the latter is studied and implemented as follows.

    CrystallineExtinction( bool has_extinction,
                           int model_option, double l, double Gg,
                           int tilt_dist_option, double L,
                           Recipe recipe,
                           bool has_texture,
                           const NCrystal::Vector& preferred_orientation1, double R1, double f1,
                           const NCrystal::Vector& preferred_orientation2, double R2, double f2,
                           const NCrystal::StructureInfo& struct_info,
                           NC::PlaneProvider * plane_provider );

    //Provide cross sections for a given neutron:
    double calcCrossSection( double neutron_ekin ) const;

    //Sample scattering event (rng is random number stream). Results are given
    //as the final ekin of the neutron and scat_mu which is cos(scattering_angle).
    struct ScatEvent { double ekin_final, mu; };
    ScatEvent sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const;

  private:
    //Data members:
    bool m_has_extinction;
    int m_model_option;
    double m_l;
    double m_Gg;
    int m_tilt_dist_option;
    double m_L;
    Recipe m_recipe;
    bool m_has_texture;
    NCrystal::Vector m_preferred_orientation1;
    double m_R1;
    double m_f1;
    NCrystal::Vector m_preferred_orientation2;
    double m_R2;
    double m_f2;
    double m_Nc;
    double m_xsectfact;
    struct HKLPlane {
      NCrystal::Vector hkl;
      double d_hkl;
      double strength;
      double F_hkl;
    };
    std::vector<HKLPlane> m_hklPlanes;
  };
  using PhysicsModel = CrystallineExtinction;

}
#endif
