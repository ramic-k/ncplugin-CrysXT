
#include "NCPluginFactory.hh"
#include "NCPhysicsModel.hh"
#include "NCrystal/internal/extd_utils/NCOrientUtils.hh"
#include "NCrystal/internal/utils/NCRotMatrix.hh"

namespace NCPluginNamespace {

  //Anisotropic (oriented) wrapper: used when the material is oriented AND has a
  //texture, so the textured coherent-elastic cross section and scattering depend on
  //the incident direction relative to the lab-fixed texture axis (MC-correct).
  class PluginScatterAniso final : public NC::ProcImpl::ScatterAnisotropicMat {
  public:
    const char * name() const noexcept override { return NCPLUGIN_NAME_CSTR "ModelAniso"; }
    PluginScatterAniso( PhysicsModel && pm ) : m_pm(std::move(pm)) {}

    NC::CrossSect crossSection( NC::CachePtr&, NC::NeutronEnergy ekin,
                                const NC::NeutronDirection& indir ) const override
    {
      return NC::CrossSect{ m_pm.calcCrossSectionDir( ekin.dbl(), indir.as<NC::Vector>().unit() ) };
    }

    NC::ScatterOutcome sampleScatter( NC::CachePtr&, NC::RNG& rng, NC::NeutronEnergy ekin,
                                      const NC::NeutronDirection& indir ) const override
    {
      auto out = m_pm.sampleScatteringEventDir( rng, ekin.dbl(), indir.as<NC::Vector>().unit() );
      return { NC::NeutronEnergy{ out.ekin_final }, out.outdir.as<NC::NeutronDirection>() };
    }
  private:
    PhysicsModel m_pm;
  };

  class PluginScatter final : public NC::ProcImpl::ScatterIsotropicMat {
  public:

    //The factory wraps our custom PhysicsModel helper class in an NCrystal API
    //Scatter class.

    const char * name() const noexcept override
    {
      return NCPLUGIN_NAME_CSTR "Model";
    }

    PluginScatter( PhysicsModel && pm ) : m_pm(std::move(pm)) {}

    //Per-neutron cache of the last (ekin -> cross section): a neutron re-queries the
    //cross section at the same energy across flights, so this avoids recomputing the
    //per-reflection sum (#5). Thread- and lifetime-safe via the NCrystal CachePtr.
    struct XSCache final : public NC::CacheBase {
      double ekin = -1.0, xs = 0.0;
      void invalidateCache() override { ekin = -1.0; }
    };

    NC::CrossSect
    crossSectionIsotropic( NC::CachePtr& cp,
                           NC::NeutronEnergy ekin ) const override
    {
      auto& c = accessCache<XSCache>(cp);
      const double e = ekin.dbl();
      if ( c.ekin != e ) { c.ekin = e; c.xs = m_pm.calcCrossSection(e); }
      return NC::CrossSect{ c.xs };
    }

    NC::ScatterOutcomeIsotropic
    sampleScatterIsotropic( NC::CachePtr&,
                            NC::RNG& rng,
                            NC::NeutronEnergy ekin ) const override
    {
      auto outcome = m_pm.sampleScatteringEvent( rng, ekin.dbl() );
      return { NC::NeutronEnergy{outcome.ekin_final},
               NC::CosineScatAngle{outcome.mu} };
    }

  private:
    PhysicsModel m_pm;
  };

}

const char * NCP::PluginFactory::name() const noexcept
{
  //Factory name. Keep this standardised form please:
  return NCPLUGIN_NAME_CSTR "Factory";
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Here follows the factory logic, for how the physics model provided by the  //
// plugin should be combined with existing models in NCrystal.                //
//                                                                            //
// In the silly example here, we want our custom physics model to replace the //
// existing incoherent-elastic model of NCrystal with our own model.          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

NC::Priority
NCP::PluginFactory::query( const NC::FactImpl::ScatterRequest& cfg ) const
{
  //Must return value >0 if we should do something, and a value higher than
  //100 means that we take precedence over the standard NCrystal factory:
  if (!cfg.get_coh_elas())
    return NC::Priority::Unable;//coherent-elastic disabled, do nothing.

  //Ok, we might be applicable. Load input data and check if is something we
  //want to handle:
  if ( ! PhysicsModel::isApplicable( cfg.info() ) )
    return NC::Priority::Unable;

  //Ok, all good. Tell the framework that we want to deal with this, with a
  //higher priority than the standard factory gives (which is 100):
  return NC::Priority{999};
}

NC::ProcImpl::ProcPtr
NCP::PluginFactory::produce( const NC::FactImpl::ScatterRequest& cfg ) const
{
  //Ok, we are selected as the provider! First create our own scatter model:

  auto sc_pp = createStdPlaneProvider( cfg.infoPtr() );
  auto pm = PhysicsModel::createFromInfo( cfg.info(), sc_pp.get() );

  //Now we just need to combine this with all the other physics.  So ask the
  //framework to set this up, except for coherent-elastic physics of course
  //since we are now dealing with that ourselves:

  auto sc_std = globalCreateScatter( cfg.modified("coh_elas=0") );

  //If the material is ORIENTED and carries a texture, provide the MC-correct,
  //direction-dependent (anisotropic) textured coherent-elastic scattering. Set the
  //crystal->lab rotation so the texture axis is fixed in the lab frame. Otherwise
  //fall back to the isotropic (orientation-averaged) model.
  if ( cfg.isSingleCrystal() && pm.hasTexture() ) {
    const auto& si = cfg.info().getStructureInfo();
    NC::RotMatrix reci = NC::getReciprocalLatticeRot( si );
    NC::RotMatrix cry2lab = NC::getCrystal2LabRot( cfg.createSCOrientation(), reci );
    pm.setOrientation( cry2lab );
    auto sc_ourmodel = NC::makeSO<PluginScatterAniso>( std::move(pm) );
    return combineProcs( sc_std, sc_ourmodel );
  }

  auto sc_ourmodel = NC::makeSO<PluginScatter>( std::move(pm) );

  //Combine and return:
  return combineProcs( sc_std, sc_ourmodel );
}
