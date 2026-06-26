#ifndef NCPlugin_CrystallineExtinction_hh
#define NCPlugin_CrystallineExtinction_hh

#include "NCrystal/NCPluginBoilerplate.hh"//Common stuff (includes NCrystal
                                          //public API headers, sets up
                                          //namespaces and aliases)
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/utils/NCRotMatrix.hh"
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

    bool hasTexture() const { return m_has_texture; }

    //Provide cross sections for a given neutron:
    double calcCrossSection( double neutron_ekin ) const;

    //Sample scattering event (rng is random number stream). Results are given
    //as the final ekin of the neutron and scat_mu which is cos(scattering_angle).
    struct ScatEvent { double ekin_final, mu; };
    ScatEvent sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const;

    //--- Anisotropic (oriented) texture path: MC-correct direction-dependent ---
    //  Enabled when the material is oriented (setOrientation called by the factory
    //  with the crystal->lab rotation). The textured coherent-elastic cross section
    //  then depends on the incident direction relative to the lab-fixed texture
    //  axis, and the scattering azimuth on each Debye cone is sampled from the
    //  modified March-Dollase pole density (model B). See NCTextureBragg notes.
    void setOrientation( const NCrystal::RotMatrix& crystal2lab );
    bool isOriented() const { return m_oriented; }

    double calcCrossSectionDir( double neutron_ekin, const NCrystal::Vector& indir ) const;

    struct ScatEventDir { double ekin_final; NCrystal::Vector outdir; };
    ScatEventDir sampleScatteringEventDir( NC::RNG& rng, double neutron_ekin,
                                           const NCrystal::Vector& indir ) const;

  private:
    //Cone-averaged textured pole density (texture factor) of plane i for an
    //incident lab direction, at wavelength wl; and the per-plane pole density at a
    //specific scattering-vector direction Qhat (lab). Both fold the (one or two)
    //March-Dollase components with weights f1,f2.
    double textureFactorDir( std::size_t iplane, const NCrystal::Vector& indir,
                             double wl ) const;
    double poleDensityAtQ( std::size_t iplane, const NCrystal::Vector& Qhat ) const;
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
      double cosA1 = 0.0, cosA2 = 0.0;     //precomputed cos(axis,normal) per component (#4)
      std::vector<double> texP;            //tabulated f-weighted texture factor vs sin(theta) (#1)
    };
    std::vector<HKLPlane> m_hklPlanes;

    //Texture-factor table: P(sin_theta) on a uniform grid sin_theta in [0,1], built once
    //at construction so the runtime cross section is a table lookup, not an integral (#1).
    static constexpr int NSINT = 1025;
    double textureFactorTab( const HKLPlane& e, double wl ) const;

    //Non-oriented total cross-section tabulation: an edge-aware grid of g(E)=xs(E)*E
    //built once at construction, so the (isotropic) calcCrossSection is a binary-search +
    //linear interpolation instead of a per-plane Bragg sum. This is O(log Ngrid) at ALL
    //energies, eliminating the O(Nplanes) cost in the epithermal regime where every plane
    //is Bragg-active. g=xs*E is used (not xs) because it is flat-1/E-free: it varies only
    //through the smooth texture/extinction factors within each inter-edge segment and
    //steps only at the Bragg edges (which become exact grid points), so linear
    //interpolation is essentially exact. The exact per-plane path is retained for
    //scattering-event sampling (which needs the per-plane weights) and for building the
    //table itself.
    double calcCrossSectionExact( double neutron_ekin ) const;
    void buildXSTable();
    std::vector<double> m_xsTabE;   //sorted neutron_ekin grid [eV]
    std::vector<double> m_xsTabG;   //g = xs*ekin at each grid point [barn*eV]

    //Oriented-texture state (set by setOrientation):
    bool m_oriented = false;
    NCrystal::Vector m_axis1_lab, m_axis2_lab;        //texture axes in lab frame (unit)
    std::vector<NCrystal::Vector> m_normal_lab;       //per-plane reflection normal in lab (unit)

    //2-D table of the modified-March-Dollase pole density md_pole_density(cosBetaQ,cosAlpha)
    //per component, precomputed at construction so the oriented (anisotropic) path is a
    //bilinear lookup instead of a per-evaluation chi-integral (#6).
    static constexpr int NBQ = 257, NAL = 129;
    std::vector<double> m_poleTab1, m_poleTab2;
    double poleTabInterp( const std::vector<double>& tab, double cosBetaQ, double cosAlpha ) const;

    //Oriented cross-section tabulation: because the texture is azimuthally symmetric about
    //each axis, the oriented cross section depends ONLY on the angle gamma_i between the
    //incident direction and texture axis i, so sigma(dir,E)=f1*sigma1(cos g1,E)+f2*sigma2(cos
    //g2,E) with each sigma_i a 2-D function. At setOrientation we tabulate g_i=sigma_i*E on a
    //(cos gamma, E) grid (E grid reused from m_xsTabE, edge-aware), making calcCrossSectionDir
    //two bilinear lookups instead of a per-plane Debye-cone average -> O(1), ~as fast as the
    //isotropic path. The exact per-plane path is kept for sampling and for building the table.
    static constexpr int NGAM = 49;                   //cos(gamma) grid points in [-1,1]
    std::vector<double> m_dirTab1, m_dirTab2;         //[NGAM * m_xsTabE.size()] of g=sigma_i*E
    std::vector<double> m_nDotA1, m_nDotA2;           //per-plane |normal . axis_i| (oriented)
    void buildOrientedXSTables();
    double extinctionFactor( const HKLPlane& e, double wl ) const;
    double calcCrossSectionDirExact( double neutron_ekin, const NCrystal::Vector& indir ) const;
    double dirTabInterp( const std::vector<double>& tab, double cosGamma, double ekin ) const;

    //3-D cone-average table G_i(cos gamma, sin theta, cos alpha) per component. The
    //Debye-cone average of the pole density depends ONLY on these three scalars
    //(cos gamma = dir.axis, sin theta = wl/2d, cos alpha = |normal.axis|), so tabulating it
    //once turns the per-plane directional texture factor into a trilinear lookup. This makes
    //both the oriented-xs table build (above) and the oriented sampler's plane selection cheap.
    static constexpr int NCG = 33, NCT = 65, NCA = 49;
    std::vector<double> m_coneTab1, m_coneTab2;       //[NCG*NCT*NCA] cone-avg pole density
    void buildConeTables();
    double coneAvgRaw( double cosGamma, double sinT, double cosAlpha,
                       const std::vector<double>& tab, unsigned npsi ) const;
    double coneTabInterp( const std::vector<double>& tab, double cosGamma, double sinT, double cosAlpha ) const;
    double textureFactorDirFast( std::size_t iplane, const NCrystal::Vector& indir, double wl ) const;
  };
  using PhysicsModel = CrystallineExtinction;

}
#endif
