////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  CrysXT: Combined crystalline extinction and texture plugin for NCrystal.  //
//                                                                            //
//  Extinction models:                                                        //
//    Sabine uncorrelated block model (Sabine_uncorr)                         //
//    Sabine correlated block model   (Sabine_corr)                           //
//    Becker & Coppens pure           (BC_pure)                               //
//    Becker & Coppens mixed          (BC_mix)                                //
//    Becker & Coppens modified       (BC_mod)                                //
//                                                                            //
//  Refs: T.M. Sabine, International Tables for Crystallography (2006),       //
//        Vol. C, Chapter 6.4, pp. 609-616.                                   //
//        P. Becker & P. Coppens, Acta Cryst. (1974). A30, 129.              //
//        P. Becker & P. Coppens, Acta Cryst. (1995). A51, 662-667.          //
//        T. Kittelmann et al., Acta Cryst. (2026). A82, 163-178.            //
//        https://doi.org/10.1107/S2053273326001245                           //
//                                                                            //
//  Texture model:                                                            //
//    Modified March-Dollase preferred orientation distribution function       //
//    with two preferred orientations and fractions f1, f2 (f1+f2=1).        //
//                                                                            //
//  Ref: H. Sato et al., J. Appl. Cryst. (2011). 44, 1128-1135.             //
//                                                                            //
//  Based on:                                                                 //
//    ncplugin-CrysExtn: https://github.com/XuShuqi7/ncplugin-CrysExtn       //
//    ncplugin-CrysText: https://github.com/highness-eu/ncplugin-CrysText     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "NCPhysicsModel.hh"
#include "NCbc2025.hh"

//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/core/NCDefs.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/extd_utils/NCOrientUtils.hh"
#include "NCrystal/internal/extd_utils/NCPlaneProvider.hh"

namespace NCPluginNamespace {

  namespace {

    // ************************************************************* //
    // ******************  Sabine's model  ************************* //
    // ************************************************************* //

    //Ref: International Tables for Crystallography (2006). Vol. C, Chapter 6.4, pp. 609–616.

    //A, B
    struct ABpair {

      double A;
      double B;

    };

    ABpair calc_AB( double y ) {

      nc_assert( y >= 0. );

      ABpair result;
      if ( y <= 1.e-9 ) {
        result.A = 1.;
        result.B = 1.;
      }
      else {
        result.A = NCrystal::exp_negarg_approx(-y) * std::sinh(y) / y;
        result.B = 1. / y - NCrystal::exp_negarg_approx(-y) / std::sinh(y);
      }

      return result;
    }

    //extinction factors
    struct ExtinctionFactors {

      double E_L; //corresponding to 2theta=0
      double E_B; //corresponding to 2theta=pi

    };

    //primary extinction factors
    ExtinctionFactors prim_extn_fact( double x, double y ) {

      ExtinctionFactors result;
      double EL = NCrystal::exp_negarg_approx(-y);
      if ( x <= 1. ) {
        EL *= (1. - x / 2. + x * x / 4. - 5. * x * x * x / 48. + 7. * x * x * x * x / 192.);
      }
      else {
        EL *= std::sqrt(2. * NCrystal::kInvPi / x );
        EL *= (1. - 1. / 8. / x - 3. / 128. / x / x - 15. / 1024. / x / x / x);
      }

      ABpair AB = calc_AB( y );

      result.E_L = EL;
      result.E_B = AB.A / std::sqrt(1. + AB.B * x);

      return result;
    }

    //secondary extinction factors
    ExtinctionFactors scnd_extn_fact( double x, double y, int tilt_dist ) {

      //tilt_dist : distribution type for the tilts between mosaic blocks
      //0 represents rectangular function, 1 for triangular function
      nc_assert( tilt_dist==0 || tilt_dist==1 );

      ExtinctionFactors result;
      ABpair AB = calc_AB( y );
      double Bx = AB.B * x;
      if ( tilt_dist == 0 ) {
        if ( x < 1.e-9 ) {
          result.E_L = NCrystal::exp_negarg_approx(-y);
        }
        else {
          result.E_L = NCrystal::exp_negarg_approx(-y) / 2. / x * (1. - NCrystal::exp_negarg_approx(-2 * x));
        }
        result.E_B = AB.A / (1. + Bx);
      }
      else {
        if ( x < 1.e-9 ) {
          result.E_L = NCrystal::exp_negarg_approx(-y);
          result.E_B = AB.A * AB.B;
        }
        else {
          result.E_L = NCrystal::exp_negarg_approx(-y) / x * (1. - (1. - NCrystal::exp_negarg_approx(-2 * x)) / 2. / x);
          result.E_B = 2 * AB.A / Bx / x * (Bx - std::log1p(Bx));
        }
      }

      return result;
    }

    //uncorrelated block model for calculating the extinction factor
    double uncorr_blk_mdl( double Nc, double wl, double F_hkl, double l,
                           double d_hkl, double mu, double G, double L, int tilt_dist ) {

      //Calculation of the extinction factor E_hkl using the uncorrelated block model
      //Nc : number of unit cells per unit volume, Aa^-3
      //wl : wavelength, Aa
      //F_hkl : |F_hkl|, modulus of the structure factor per unit cell, Aa
      //l : block size, Aa
      //d_hkl : dspacing for the hkl plan, Aa
      //mu : incoherent and absorption cross section per unit volume, Aa^-1
      //G : integral breadth of the angular distribution of mosaic blocks, dimensionless
      //L : grain size, Aa (A grain is formed by crystallites or blocks.)
      //tilt_dist : distribution type for the tilts between mosaic blocks
      //0 represents rectangular function, 1 for triangular function

      double sin_theta = 0.5 * wl / d_hkl; //2*d_hkl*sin(theta_hkl)=lambda
      if ( sin_theta <= 1. ) {
        double sin_theta_square = NC::ncsquare(sin_theta);
        double cos_theta_square = 1. - sin_theta_square;
        double cos_theta = std::sqrt(cos_theta_square);
        double y = mu * l;

        //primary extinction
        double xp = NC::ncsquare(Nc * wl * F_hkl * l);
        ExtinctionFactors EpLB = prim_extn_fact( xp, y );
        double Ep = EpLB.E_L * cos_theta_square + EpLB.E_B * sin_theta_square;

        //secondary extinction
        if ( sin_theta != 0. && cos_theta != 0. ) {
          double Q_theta = NC::ncsquare(Nc * wl * F_hkl) * wl / 2. / sin_theta / cos_theta;
          //double xs = Ep * Q_theta * G * l;
          double xs = Ep * Q_theta * G * L; //from equation (6.4.9.1) in "International Tables for Crystallography (2006). Vol. C, Chapter 6.4, pp. 609–616."
          ExtinctionFactors EsLB = scnd_extn_fact( xs, y, tilt_dist );
          double Es = EsLB.E_L * cos_theta_square + EsLB.E_B * sin_theta_square;

          return Ep * Es;
        }
        else {
          return 0.;
        }
      }
      else {

        return 1.;
      }
    }

    //correlated block model for calculating the extinction factor
    double corr_blk_mdl( double Nc, double wl, double F_hkl, double l,
                         double d_hkl, double mu, double g, double L ) {

      //Calculation of the extinction factor E_hkl using the uncorrelated block model
      //Nc : number of unit cells per unit volume, Aa^-3
      //wl : wavelength, Aa
      //F_hkl : |F_hkl|, modulus of the structure factor per unit cell, Aa
      //l : block size, Aa
      //d_hkl : dspacing for the hkl plan, Aa
      //mu : incoherent and absorption cross section per unit volume, Aa^-1
      //g : standard deviation of the distribution of tilts * \sqrt(pi)/2, dimensionless
      //L : side of cube of the crystal, Aa (meaning to be clarified)

      double sin_theta = 0.5 * wl / d_hkl; //2*d_hkl*sin(theta_hkl)=lambda
      if ( sin_theta <= 1. ) {
        double sin_theta_square = NC::ncsquare(sin_theta);
        double cos_theta_square = 1. - sin_theta_square;
        double cos_theta = std::sqrt(cos_theta_square);
        double y = mu * l;

        //refine both primary and secondary extinction in this model
        double E;
        if ( l > 0. && g == 0. ) {
          double x = NC::ncsquare(Nc * wl * F_hkl * l); //pure primary
          ExtinctionFactors ELB = prim_extn_fact( x, y );
          E = ELB.E_L * cos_theta_square + ELB.E_B * sin_theta_square;
        }
        else {
          if ( sin_theta != 0. && cos_theta != 0. ) {
            double Q_theta = NC::ncsquare(Nc * wl * F_hkl) * wl / 2. / sin_theta / cos_theta;
            double x = NC::ncsquare(Nc * wl * F_hkl * l + g * Q_theta * (L - l));
            ExtinctionFactors ELB = prim_extn_fact( x, y );
            E = ELB.E_L * cos_theta_square + ELB.E_B * sin_theta_square;
          }
          else {
            E = 0.;
          }
        }

        return E;
      }
      else {

        return 1.;
      }
    }

    // ************************************************************* //
    // ***************  Becker & Coppens' model ******************** //
    // ************************************************************* //

    // Refs: Acta Cryst. (1974). A30, 129
    //       Acta Cryst. (1995). A51, 662-667

    ABpair calc_AB_theta( double cos_2theta, int opt ) {

      //Calculations of A(theta) and B(theta)
      //cos_2theta : cos(2*theta_hkl)
      //opt : 0 for primary extinction, 1, 2, 3 for sencondary extinction following a
      //Gaussian, Lorentzian or Fresnel distribution, respectively

      nc_assert( opt <= 3 );

      ABpair result;

      if ( opt == 0 ) {
        result.A = 0.20 + 0.45 * cos_2theta;
        result.B = 0.22 - 0.12 * NC::ncsquare(0.5 - cos_2theta);
      }
      else if ( opt == 1 ) {
        result.A = 0.58 + 0.48 * cos_2theta + 0.24 * NC::ncsquare(cos_2theta);
        result.B = 0.02 - 0.025 * cos_2theta;
      }
      else if ( opt == 2 ) {
        result.A = 0.025 + 0.285 * cos_2theta;
        if ( cos_2theta >= 0. ) {
          result.B = 0.15 - 0.2 * NC::ncsquare(0.75 - cos_2theta);
        }
        else {
          result.B = -0.45 * cos_2theta;
        }
      }
      else {
        result.A = 0.48 + 0.6 * cos_2theta;
        result.B = 0.20 - 0.06 * NC::ncsquare(0.2 - cos_2theta);
      }

      return result;
    }

    double BC_pure_extn_mdl( double Nc, double wl, double F_hkl, double l,
                             double d_hkl, double g, double L, int tilt_dist, Recipe recipe ) {

      //Calculation of pure primary or secondary extinction factor y using the model of Becker & Coppens
      // pure primary:  BC_pure  l
      // pure secondary type-I:  BC_pure  g  L
      // pure secondary type-II:  BC_pure  l  L  Gauss/Lorentz/Fresnel
      //Nc : number of unit cells per unit volume, Aa^-3
      //wl : wavelength, Aa
      //F_hkl : |F_hkl|, modulus of the structure factor per unit cell, Aa
      //l : "t", mean path length through a perfect crystal, equivalent to block size, Aa
      //d_hkl : dspacing for the hkl plan, Aa
      //g : width parameter of the mosaic distribution, dimensionless
      //L : "T-bar", mean path length through a mosaic crystal, Aa
      //tilt_dist : option for the calculation of A(theta) and B(theta),
      //1, 2, 3 for orientation of crystallite following a Gaussian, Lorentzian
      //or Fresnel distribution, respectively

      double sin_theta = 0.5 * wl / d_hkl; //2*d_hkl*sin(theta_hkl)=lambda
      if ( sin_theta <= 1. ) {
        double cos_theta  = std::sqrt(1. - NC::ncsquare(sin_theta));
        double sin_2theta = 2. * sin_theta * cos_theta;
        double cos_2theta = 1. - 2. * NC::ncsquare(sin_theta);
        //double Q_theta = NC::ncsquare(Nc * wl * F_hkl) * wl / sin_2theta; //same as in Sabine's model
        double Q_theta = NC::ncsquare(Nc * wl * F_hkl) * wl; //division by sin_2theta to be done later

        double y;
        //pure primary extinction
        if ( l > 0. && g == 0. && L == 0. ) {
          //double x = 2. / 3. * Q_theta * l * l * sin_2theta / wl;
          double x = 2. / 3. * Q_theta * l * l / wl;
          if ( recipe == Recipe::cls ) {
            ABpair AB_theta = calc_AB_theta( cos_2theta, 0 );
            y = 1. / std::sqrt(1. + 2. * x + AB_theta.A * NC::ncsquare(x) / (1. + AB_theta.B * x));
          } else if ( recipe == Recipe::std ) {
            y = bc2025_y_primary( x, sin_theta );
          } else {
            y = bc2025_y_primary_lux( x, sin_theta );
          }
        }
        //pure secondary extinction type-I
        else if ( l == 0. && g > 0. && L > 0. ) {
          //double x = 2. * std::sqrt(2.) / 3. * g * Q_theta * L;
          if ( sin_2theta != 0. ) {
            Q_theta /= sin_2theta;
            y = 1. / std::sqrt(1. + 2. * g * Q_theta * L);
          }
          else {
            y = 0.;
          }
        }
        //pure secondary extinction type-II
        else if ( l > 0. && g == 0. && L > 0. ) {
          //double x = 2. / 3. * Q_theta * L * l * sin_2theta / wl;
          double x = 2. / 3. * Q_theta * L * l / wl;
          if ( recipe == Recipe::cls ) {
            ABpair AB_theta = calc_AB_theta( cos_2theta, tilt_dist );
            y = 1. / std::sqrt(1. + 2.12 * x + AB_theta.A * NC::ncsquare(x) / (1. + AB_theta.B * x));
          } else if ( recipe == Recipe::std ) {
            if      ( tilt_dist == 1 ) y = bc2025_y_scndgauss(   x, sin_theta );
            else if ( tilt_dist == 2 ) y = bc2025_y_scndlorentz(  x, sin_theta );
            else                       y = bc2025_y_scndfresnel(   x, sin_theta );
          } else {
            if      ( tilt_dist == 1 ) y = bc2025_y_scndgauss_lux(   x, sin_theta );
            else if ( tilt_dist == 2 ) y = bc2025_y_scndlorentz_lux(  x, sin_theta );
            else                       y = bc2025_y_scndfresnel_lux(   x, sin_theta );
          }
        }
        else {
          y = 1.;
        }

        return y;
      }
      else {

        return 1.;
      }

    }

    double BC_mix_extn_mdl( double Nc, double wl, double F_hkl, double l,
                            double d_hkl, double g, double L, int tilt_dist, Recipe recipe ) {

      //Calculation of mixed primary or secondary extinction factor y using the model of Becker & Coppens
      //Nc : number of unit cells per unit volume, Aa^-3
      //wl : wavelength, Aa
      //F_hkl : |F_hkl|, modulus of the structure factor per unit cell, Aa
      //l : "t", mean path length through a perfect crystal, equivalent to block size, Aa
      //d_hkl : dspacing for the hkl plan, Aa
      //g : width parameter of the mosaic distribution, dimensionless
      //L : "T-bar", mean path length through a mosaic crystal, Aa
      //tilt_dist : option for the calculation of A(theta) and B(theta),
      //1, 2, 3 for orientation of crystallite following a Gaussian, Lorentzian
      //or Fresnel distribution, respectively

      double sin_theta = 0.5 * wl / d_hkl; //2*d_hkl*sin(theta_hkl)=lambda
      if ( sin_theta >= -1. && sin_theta <= 1. ) {
        double cos_theta  = std::sqrt(1. - NC::ncsquare(sin_theta));
        double sin_2theta = 2. * sin_theta * cos_theta;
        double cos_2theta = 1. - 2. * NC::ncsquare(sin_theta);
        //double Q_theta = NC::ncsquare(Nc * wl * F_hkl) * wl / sin_2theta; //same as in Sabine's model
        double Q_theta = NC::ncsquare(Nc * wl * F_hkl) * wl; //division by sin_2theta to be done later

        //primary
        //double xp = 2. / 3. * Q_theta * l * l * sin_2theta / wl;
        double xp = 2. / 3. * Q_theta * l * l / wl;
        double yp;
        if ( recipe == Recipe::cls ) {
          ABpair AB_theta_p = calc_AB_theta( cos_2theta, 0 );
          yp = 1. / std::sqrt(1. + 2. * xp + AB_theta_p.A * NC::ncsquare(xp) / (1. + AB_theta_p.B * xp));
        } else if ( recipe == Recipe::std ) {
          yp = bc2025_y_primary( xp, sin_theta );
        } else {
          yp = bc2025_y_primary_lux( xp, sin_theta );
        }

        double xs, ys;
        if ( l < 1.e-9 ) {
          xs = 0.;
          ys = 1.; //becomes pure primary
        }
        else {
          //xs = 2. / 3. * Q_theta * L / std::sqrt(NC::ncsquare(wl / l / sin_2theta) + 1. / (2. * g * g));
		  if ( tilt_dist == 1 || tilt_dist == 3 ) {
			xs = 2. / 3. * Q_theta * L / std::sqrt(NC::ncsquare(wl / l) + NC::ncsquare(sin_2theta) / (2. * g * g)); //Valid for Gaussian and Fresnel distributions, Eq.40(b) in Acta Cryst. (1974). A30, 129
		  }
		  else {
            xs = 2. / 3. * Q_theta * L / (wl / l + sin_2theta * 2. / (3. * g)); //Valid for Lorentzian distribution, Eq.41(b) in Acta Cryst. (1974). A30, 129
		  }
          xs *= yp; //Correction of formula, ys also dependent on yp
          if ( recipe == Recipe::cls ) {
            ABpair AB_theta_s = calc_AB_theta( cos_2theta, tilt_dist );
		    if ( tilt_dist == 1 ) {
			  ys = 1. / std::sqrt(1. + 2.12 * xs + AB_theta_s.A * NC::ncsquare(xs) / (1. + AB_theta_s.B * xs)); //The factor 2.12 is only applied in the case of Gaussian distribution
		    }
		    else {
              ys = 1. / std::sqrt(1. + 2 * xs + AB_theta_s.A * NC::ncsquare(xs) / (1. + AB_theta_s.B * xs));
		    }
          } else if ( recipe == Recipe::std ) {
            if      ( tilt_dist == 1 ) ys = bc2025_y_scndgauss(   xs, sin_theta );
            else if ( tilt_dist == 2 ) ys = bc2025_y_scndlorentz(  xs, sin_theta );
            else                       ys = bc2025_y_scndfresnel(   xs, sin_theta );
          } else {
            if      ( tilt_dist == 1 ) ys = bc2025_y_scndgauss_lux(   xs, sin_theta );
            else if ( tilt_dist == 2 ) ys = bc2025_y_scndlorentz_lux(  xs, sin_theta );
            else                       ys = bc2025_y_scndfresnel_lux(   xs, sin_theta );
          }
        }

        return yp * ys;
      }
      else {

        return 1.;
      }
    }

    double BC_mod_extn_mdl( double Nc, double wl, double F_hkl, double l,
                            double d_hkl, double g, double L, int tilt_dist, Recipe recipe ) {

      //Calculation of mixed primary or secondary extinction factor y using the MODIFIED model of Becker & Coppens
      //Only SECONDARY extinction can happen, but is characterized by l and g
      //Nc : number of unit cells per unit volume, Aa^-3
      //wl : wavelength, Aa
      //F_hkl : |F_hkl|, modulus of the structure factor per unit cell, Aa
      //l : "t", mean path length through a perfect crystal, equivalent to block size, Aa
      //d_hkl : dspacing for the hkl plan, Aa
      //g : width parameter of the mosaic distribution, dimensionless
      //L : "T-bar", mean path length through a mosaic crystal, Aa
      //tilt_dist : option for the calculation of A(theta) and B(theta),
      //1, 2, 3 for orientation of crystallite following a Gaussian, Lorentzian
      //or Fresnel distribution, respectively

      double sin_theta = 0.5 * wl / d_hkl; //2*d_hkl*sin(theta_hkl)=lambda
      if ( sin_theta >= -1. && sin_theta <= 1. ) {
        double cos_theta  = std::sqrt(1. - NC::ncsquare(sin_theta));
        double sin_2theta = 2. * sin_theta * cos_theta;
        double cos_2theta = 1. - 2. * NC::ncsquare(sin_theta);
        //double Q_theta = NC::ncsquare(Nc * wl * F_hkl) * wl / sin_2theta; //same as in Sabine's model
        double Q_theta = NC::ncsquare(Nc * wl * F_hkl) * wl; //division by sin_2theta to be done later

        //primary
        //double xp = 2. / 3. * Q_theta * l * l * sin_2theta / wl;
        //double xp = 2. / 3. * Q_theta * l * l / wl;
        //double yp = 1. / std::sqrt(1. + 2. * xp + AB_theta_p.A * NC::ncsquare(xp) / (1. + AB_theta_p.B * xp));
        double yp = 1.; //No primary extinction

        double xs, ys;
        if ( l < 1.e-9 ) {
          xs = 0.;
          ys = 1.; //becomes pure primary
        }
        else {
          //xs = 2. / 3. * Q_theta * L / std::sqrt(NC::ncsquare(wl / l / sin_2theta) + 1. / (2. * g * g));
		  if ( tilt_dist == 1 || tilt_dist == 3 ) {
			xs = 2. / 3. * Q_theta * L / std::sqrt(NC::ncsquare(wl / l) + NC::ncsquare(sin_2theta) / (2. * g * g)); //Valid for Gaussian and Fresnel distributions, Eq.40(b) in Acta Cryst. (1974). A30, 129
		  }
		  else {
            xs = 2. / 3. * Q_theta * L / (wl / l + sin_2theta * 2. / (3. * g)); //Valid for Lorentzian distribution, Eq.41(b) in Acta Cryst. (1974). A30, 129
		  }
          xs *= yp; //Correction of formula, ys also dependent on yp
          if ( recipe == Recipe::cls ) {
            ABpair AB_theta_s = calc_AB_theta( cos_2theta, tilt_dist );
		    if ( tilt_dist == 1 ) {
			  ys = 1. / std::sqrt(1. + 2.12 * xs + AB_theta_s.A * NC::ncsquare(xs) / (1. + AB_theta_s.B * xs)); //The factor 2.12 is only applied in the case of Gaussian distribution
		    }
		    else {
              ys = 1. / std::sqrt(1. + 2 * xs + AB_theta_s.A * NC::ncsquare(xs) / (1. + AB_theta_s.B * xs));
		    }
          } else if ( recipe == Recipe::std ) {
            if      ( tilt_dist == 1 ) ys = bc2025_y_scndgauss(   xs, sin_theta );
            else if ( tilt_dist == 2 ) ys = bc2025_y_scndlorentz(  xs, sin_theta );
            else                       ys = bc2025_y_scndfresnel(   xs, sin_theta );
          } else {
            if      ( tilt_dist == 1 ) ys = bc2025_y_scndgauss_lux(   xs, sin_theta );
            else if ( tilt_dist == 2 ) ys = bc2025_y_scndlorentz_lux(  xs, sin_theta );
            else                       ys = bc2025_y_scndfresnel_lux(   xs, sin_theta );
          }
        }

        return yp * ys;
      }
      else {

        return 1.;
      }
    }

    // ************************************************************* //
    // ******************  Texture model  ************************* //
    // ************************************************************* //

    //Ref: Sato et al. 2011

    double cosAngleVectors( const NC::Vector& a, const NC::Vector& b )
    {
      return NC::ncclamp( a.dot(b) / ( std::sqrt( a.mag2() * b.mag2() ) ),
                          -1.0, 1.0 );
    }

    //preferred orientation distribution function (Sato 2011)
    double sato_mmd_podf( const NCrystal::Vector& preferred_orientation, NCrystal::Vector vec_hkl,
                          double d_hkl, double R, double wl )
    {
      double P_hkl = 1.;
      const unsigned int num_phis = 1000;

      double sin_theta = 0.5 * wl / d_hkl;
      if ( sin_theta >= -1. && sin_theta <= 1. ) {
        double cos_theta = std::sqrt( 1.0 - NC::ncsquare(sin_theta) );

        double cos_A = cosAngleVectors( preferred_orientation, vec_hkl );
        double sin_A = std::sqrt( 1.0 - NC::ncsquare(cos_A) );

        //trapezoidal integration
        P_hkl = 0.;
        for ( auto phi : NC::linspace( 0, NC::k2Pi * (1-1./num_phis), num_phis ) ) {
          double B = cos_A * sin_theta + sin_A * cos_theta * std::sin(phi);
          P_hkl += std::pow( (NC::ncsquare(R * B) + (1.0 - NC::ncsquare(B)) / R), -1.5 ) / (num_phis+1);
        }
      }

      return P_hkl;
    }

  }
}
// ************************************************************* //
// *****************  parsing & processing  ******************** //
// ************************************************************* //

bool NCP::CrystallineExtinction::isApplicable( const NC::Info& info ) {

  //Accept if input is NCMAT data with @CUSTOM_<pluginname> section:
  return info.countCustomSections(pluginNameUpperCase()) > 0;
}

NCP::CrystallineExtinction NCP::CrystallineExtinction::createFromInfo( const NC::Info& info,
                                                                       NC::PlaneProvider * plane_provider )
{
  if ( info.countCustomSections( pluginNameUpperCase() ) != 1 )
    NCRYSTAL_THROW2(BadInput,"Multiple @CUSTOM_"<<pluginNameUpperCase()<<" sections are not allowed");
  auto data = info.getCustomSection( pluginNameUpperCase() );

  // @CUSTOM_CRYSXT accepts lines beginning with the keywords Extinction or Texture:
  //
  //   Extinction  <model>  l  g  L  <dist>
  //   Texture     px py pz  R1  f1
  //   Texture     px py pz  R2  f2
  //
  // At least one keyword must be present. Texture requires exactly two lines.

  bool has_extinction = false;
  bool has_texture    = false;

  int model_option     = 0;
  double l = 0., Gg = 0., L = 0.;
  int tilt_dist_option = 0;
  Recipe recipe        = Recipe::std;

  NCrystal::Vector preferred_orientation1, preferred_orientation2;
  double R1 = 1., f1 = 1., R2 = 1., f2 = 0.;

  int texture_lines_found = 0;

  for ( auto& line : data ) {
    if ( line.empty() )
      continue;

    if ( line.at(0).compare("Extinction") == 0 ) {

      if ( has_extinction )
        NCRYSTAL_THROW2(BadInput,"Multiple Extinction lines in @CUSTOM_"<<pluginNameUpperCase()<<" are not allowed");
      has_extinction = true;

      if ( line.size() < 5 || line.size() > 7 )
        NCRYSTAL_THROW2(BadInput,"Extinction line in @CUSTOM_"<<pluginNameUpperCase()
                        <<" should have 5-7 entries (Extinction model l g L [dist] [rec=cls|std|lux])");

      if ( line.at(1).compare("Sabine_uncorr") == 0 ) {
        model_option = 0;
      }
      else if ( line.at(1).compare("Sabine_corr") == 0 ) {
        model_option = 1;
      }
      else if ( line.at(1).compare("BC_pure") == 0 ) {
        model_option = 2;
      }
      else if ( line.at(1).compare("BC_mix") == 0 ) {
        model_option = 3;
      }
      else if ( line.at(1).compare("BC_mod") == 0 ) {
        model_option = 6;
      }
      else {
        NCRYSTAL_THROW2(BadInput,"Only Sabine_uncorr, Sabine_corr, BC_pure, BC_mix and BC_mod are supported.");
      }

      if (   ! NC::safe_str2dbl( line.at(2), l  )
             || ! NC::safe_str2dbl( line.at(3), Gg )
             || ! NC::safe_str2dbl( line.at(4), L  )
             || ! (l  >= 0.0)
             || ! (Gg >= 0.0)
             || ! (L  >= 0.0) )
        NCRYSTAL_THROW2( BadInput,"Invalid values in Extinction line: l, g and L should be non-negative." );

      if ( model_option == 0 && line.size() != 6 ) {
        NCRYSTAL_THROW2(BadInput,"Extinction line for Sabine_uncorr requires exactly 6 entries (Extinction model l g L dist).");
      } else if ( (model_option == 2 || model_option == 3 || model_option == 6) && line.size() != 6 && line.size() != 7 ) {
        NCRYSTAL_THROW2(BadInput,"Extinction line for BC models requires 6 entries (with dist) or 7 entries (with dist and rec=cls|std|lux).");
      }

      if ( model_option == 0 ) {
        if ( line.at(5).compare("rect") == 0 ) {
          tilt_dist_option = 0;
        }
        else if ( line.at(5).compare("tri") == 0 ) {
          tilt_dist_option = 1;
        }
        else {
          NCRYSTAL_THROW2( BadInput,"Distribution option for Sabine_uncorr should be rect or tri." );
        }
      }
      else if ( model_option == 2 || model_option == 3 || model_option == 6 ) {
        if ( line.at(5).compare("Gauss") == 0 ) {
          tilt_dist_option = 1;
        }
        else if ( line.at(5).compare("Lorentz") == 0 ) {
          tilt_dist_option = 2;
        }
        else if ( line.at(5).compare("Fresnel") == 0 ) {
          tilt_dist_option = 3;
        }
        else {
          NCRYSTAL_THROW2( BadInput,"Distribution option for BC models should be Gauss, Lorentz or Fresnel." );
        }
        if ( line.size() == 7 ) {
          if ( line.at(6).compare("rec=cls") == 0 ) {
            recipe = Recipe::cls;
          }
          else if ( line.at(6).compare("rec=std") == 0 ) {
            recipe = Recipe::std;
          }
          else if ( line.at(6).compare("rec=lux") == 0 ) {
            recipe = Recipe::lux;
          }
          else {
            NCRYSTAL_THROW2( BadInput,"Invalid recipe option '"<<line.at(6)
                             <<"'. Expected rec=cls, rec=std, or rec=lux." );
          }
        }
      }

    }
    else if ( line.at(0).compare("Texture") == 0 ) {

      if ( line.size() != 6 )
        NCRYSTAL_THROW2(BadInput,"Texture line in @CUSTOM_"<<pluginNameUpperCase()
                        <<" should have six entries (Texture px py pz R f)");

      ++texture_lines_found;
      if ( texture_lines_found == 1 ) {
        if ( ! NC::safe_str2dbl( line.at(1), preferred_orientation1.at(0) )
             || ! NC::safe_str2dbl( line.at(2), preferred_orientation1.at(1) )
             || ! NC::safe_str2dbl( line.at(3), preferred_orientation1.at(2) )
             || ! NC::safe_str2dbl( line.at(4), R1 )
             || ! NC::safe_str2dbl( line.at(5), f1 )
             || ! (preferred_orientation1.mag()>0) || ! (R1>0.0) || ! (f1>0.0) )
          NCRYSTAL_THROW2( BadInput,"Invalid values in first Texture line." );
      }
      else if ( texture_lines_found == 2 ) {
        if ( ! NC::safe_str2dbl( line.at(1), preferred_orientation2.at(0) )
             || ! NC::safe_str2dbl( line.at(2), preferred_orientation2.at(1) )
             || ! NC::safe_str2dbl( line.at(3), preferred_orientation2.at(2) )
             || ! NC::safe_str2dbl( line.at(4), R2 )
             || ! NC::safe_str2dbl( line.at(5), f2 )
             || ! (preferred_orientation2.mag()>0) || ! (R2>0.0) || ! (f2>0.0) )
          NCRYSTAL_THROW2( BadInput,"Invalid values in second Texture line." );
      }
      else {
        NCRYSTAL_THROW2(BadInput,"More than two Texture lines in @CUSTOM_"<<pluginNameUpperCase()<<" are not allowed");
      }

    }
    else {
      NCRYSTAL_THROW2(BadInput,"Unknown keyword '"<<line.at(0)<<"' in @CUSTOM_"<<pluginNameUpperCase()
                      <<". Expected Extinction or Texture.");
    }
  }

  if ( !has_extinction && texture_lines_found == 0 )
    NCRYSTAL_THROW2(BadInput,"@CUSTOM_"<<pluginNameUpperCase()<<" section is empty or has no valid entries.");

  if ( texture_lines_found == 1 )
    NCRYSTAL_THROW2(BadInput,"Texture requires exactly two lines in @CUSTOM_"<<pluginNameUpperCase());

  if ( texture_lines_found == 2 ) {
    has_texture = true;
    if ( !(f1+f2==1.0) )
      NCRYSTAL_THROW2( BadInput,"Texture f1 and f2 must sum to 1." );
  }

  if ( !info.hasStructureInfo() )
    NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks Structure information.");
  const NCrystal::StructureInfo& struct_info = info.getStructureInfo();

  return CrystallineExtinction( has_extinction, model_option, l, Gg, tilt_dist_option, L,
                                recipe,
                                has_texture, preferred_orientation1, R1, f1,
                                preferred_orientation2, R2, f2,
                                struct_info, plane_provider );
}

NCP::CrystallineExtinction::CrystallineExtinction( bool has_extinction,
                                                   int model_option, double l, double Gg,
                                                   int tilt_dist_option, double L,
                                                   Recipe recipe,
                                                   bool has_texture,
                                                   const NCrystal::Vector& preferred_orientation1, double R1, double f1,
                                                   const NCrystal::Vector& preferred_orientation2, double R2, double f2,
                                                   const NCrystal::StructureInfo& struct_info,
                                                   NC::PlaneProvider * plane_provider )
: m_has_extinction(has_extinction),
  m_model_option(model_option),
  m_l(l),
  m_Gg(Gg),
  m_tilt_dist_option(tilt_dist_option),
  m_L(L),
  m_recipe(recipe),
  m_has_texture(has_texture),
  m_preferred_orientation1(preferred_orientation1),
  m_R1(R1),
  m_f1(f1),
  m_preferred_orientation2(preferred_orientation2),
  m_R2(R2),
  m_f2(f2)
{
  nc_assert( m_l  >= 0.0 );
  nc_assert( m_Gg >= 0.0 );
  nc_assert( m_L  >= 0.0 );
  nc_assert( plane_provider != nullptr );
  nc_assert( plane_provider->canProvide() );

  m_Nc = 1. / struct_info.volume;
  m_xsectfact = 0.5 / struct_info.n_atoms / struct_info.volume;

  NCrystal::RotMatrix lattice_rot = NC::getLatticeRot( struct_info.lattice_a, struct_info.lattice_b, struct_info.lattice_c,
                                                       struct_info.alpha*NC::kDeg, struct_info.beta*NC::kDeg, struct_info.gamma*NC::kDeg );

  plane_provider->prepareLoop();
  NCrystal::Optional<NC::PlaneProvider::Plane> opt_plane;
  while ( ( opt_plane = plane_provider->getNextPlane() ).has_value() ) {
    auto& pl = opt_plane.value();
    nc_assert( pl.dspacing > 0.0 );
    m_hklPlanes.push_back( HKLPlane{} );
    auto& e = m_hklPlanes.back();
    e.hkl      = lattice_rot * pl.demi_normal;
    e.d_hkl    = pl.dspacing;
    e.strength = pl.dspacing * pl.fsq * m_xsectfact;
    e.F_hkl    = std::sqrt(pl.fsq) * 1.e-4;
  }
}

double NCP::CrystallineExtinction::calcCrossSection( double neutron_ekin ) const {

  double xs_in_barns = 0.0;
  const double wl    = NC::ekin2wl( neutron_ekin );
  const double wlsq  = NC::ncsquare( wl );
  const double mu    = 0.;

  for ( auto& e : m_hklPlanes ) {

    if ( wl > 2 * e.d_hkl )
      break;

    double E_hkl = 1.0;
    if ( m_has_extinction ) {
      if ( m_model_option == 0 ) {
        E_hkl = uncorr_blk_mdl( m_Nc, wl, e.F_hkl, m_l, e.d_hkl, mu, m_Gg, m_L, m_tilt_dist_option );
      }
      else if ( m_model_option == 1 ) {
        E_hkl = corr_blk_mdl( m_Nc, wl, e.F_hkl, m_l, e.d_hkl, mu, m_Gg, m_L );
      }
      else if ( m_model_option == 2 ) {
        E_hkl = BC_pure_extn_mdl( m_Nc, wl, e.F_hkl, m_l, e.d_hkl, m_Gg, m_L, m_tilt_dist_option, m_recipe );
      }
      else if ( m_model_option == 3 ) {
        E_hkl = BC_mix_extn_mdl( m_Nc, wl, e.F_hkl, m_l, e.d_hkl, m_Gg, m_L, m_tilt_dist_option, m_recipe );
      }
      else {
        E_hkl = BC_mod_extn_mdl( m_Nc, wl, e.F_hkl, m_l, e.d_hkl, m_Gg, m_L, m_tilt_dist_option, m_recipe );
      }
    }

    double texture_factor = 1.0;
    if ( m_has_texture ) {
      double P1 = sato_mmd_podf( m_preferred_orientation1, e.hkl, e.d_hkl, m_R1, wl );
      double P2 = sato_mmd_podf( m_preferred_orientation2, e.hkl, e.d_hkl, m_R2, wl );
      texture_factor = P1 * m_f1 + P2 * m_f2;
    }

    xs_in_barns += e.strength * E_hkl * texture_factor;
  }
  xs_in_barns *= 2. * wlsq;

  return xs_in_barns;
}

NCP::CrystallineExtinction::ScatEvent NCP::CrystallineExtinction::sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const {

  // NOTE: This sampling does not include the effect of texture and thus this plugin should not be used for Monte-Carlo simulations if texture is present.

  ScatEvent result;
  result.ekin_final = neutron_ekin;

  const double wl   = NC::ekin2wl(neutron_ekin);
  const double wlsq = NC::ncsquare(wl);
  const double xs   = calcCrossSection( neutron_ekin ) / ( 2. * wlsq );
  const double mu   = 0.;
  const double rnd  = rng.generate();

  double left_bound  = 0.;
  double right_bound = 0.;

  for ( auto& e : m_hklPlanes ) {

    if ( wl > 2 * e.d_hkl )
      break;

    double E_hkl = 1.0;
    if ( m_has_extinction ) {
      if ( m_model_option == 0 ) {
        E_hkl = uncorr_blk_mdl( m_Nc, wl, e.F_hkl, m_l, e.d_hkl, mu, m_Gg, m_L, m_tilt_dist_option );
      }
      else if ( m_model_option == 1 ) {
        E_hkl = corr_blk_mdl( m_Nc, wl, e.F_hkl, m_l, e.d_hkl, mu, m_Gg, m_L );
      }
      else if ( m_model_option == 2 ) {
        E_hkl = BC_pure_extn_mdl( m_Nc, wl, e.F_hkl, m_l, e.d_hkl, m_Gg, m_L, m_tilt_dist_option, m_recipe );
      }
      else if ( m_model_option == 3 ) {
        E_hkl = BC_mix_extn_mdl( m_Nc, wl, e.F_hkl, m_l, e.d_hkl, m_Gg, m_L, m_tilt_dist_option, m_recipe );
      }
      else {
        E_hkl = BC_mod_extn_mdl( m_Nc, wl, e.F_hkl, m_l, e.d_hkl, m_Gg, m_L, m_tilt_dist_option, m_recipe );
      }
    }

    right_bound += e.strength * E_hkl / xs;
    nc_assert( left_bound < right_bound && right_bound <= 1.0 );

    if ( left_bound <= rnd && right_bound > rnd ) {
      const double En_hkl = 0.5 * NC::kPiSq * NC::const_hhm / NC::ncsquare(e.d_hkl);
      const double mu_n   = 1. - 2 * En_hkl / neutron_ekin;
      nc_assert( NC::ncabs(mu_n) <= 1.0 );
      result.mu = mu_n;
      break;
    }
    else {
      left_bound = right_bound;
    }
  }

  return result;
}
