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
#include <algorithm>
#include <cmath>
#include <vector>

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

    //Sato 2011 modified March-Dollase pole factor for one reflection, as a function
    //of the (signed) cosine cos_A between the texture axis and the reflection normal
    //and the Bragg sine sin_theta = wl/(2d). The azimuthal integrand is smooth and
    //2*pi-periodic, so the midpoint rule converges EXPONENTIALLY: 64 points already
    //reach ~machine precision (this is why 1000 was wasteful). pow(x,-1.5) is replaced
    //by 1/(x*sqrt(x)). This kernel is only evaluated at construction (see the per-plane
    //tabulation in the CrystallineExtinction constructor); runtime is a table lookup.
    double sato_kernel( double cos_A, double sin_theta, double R, unsigned num_phis = 64 )
    {
      if ( !( sin_theta >= -1. && sin_theta <= 1. ) )
        return 1.0;
      const double cos_theta = std::sqrt( std::max(0.0, 1.0 - sin_theta*sin_theta) );
      const double sin_A = std::sqrt( std::max(0.0, 1.0 - cos_A*cos_A) );
      double P = 0.;
      for ( unsigned i=0; i<num_phis; ++i ) {
        const double phi = NC::k2Pi * ( i + 0.5 ) / num_phis;
        const double B = cos_A * sin_theta + sin_A * cos_theta * std::sin(phi);
        const double x = R*R*B*B + ( 1.0 - B*B ) / R;
        P += 1.0 / ( x * std::sqrt(x) );
      }
      return P / num_phis;
    }

    // ********  Anisotropic (oriented) textured coherent-elastic: model B  ******** //
    //
    // Modified March-Dollase POLE DENSITY of a reflection whose normal is at angle
    // alpha to a texture axis, evaluated for a scattering vector at angle betaQ to
    // that axis: the average of the fundamental March-Dollase P_MD over the
    // crystallite free rotation chi about the scattering vector. Reduces to P_MD for
    // alpha=0, and to 1 for R=1. (Validated to machine precision against sato_mmd_podf
    // for betaQ=pi/2-theta in the Python prototype.)
    double md_pole_density( double cos_betaQ, double cos_alpha, double R, unsigned nchi = 64 )
    {
      const double cb = NC::ncclamp(cos_betaQ,-1.0,1.0);
      const double ca = NC::ncclamp(cos_alpha,-1.0,1.0);
      const double sb = std::sqrt(std::max(0.0,1.0-cb*cb));
      const double sa = std::sqrt(std::max(0.0,1.0-ca*ca));
      double sum = 0.0;
      for ( unsigned i=0; i<nchi; ++i ) {
        const double chi = NC::k2Pi * ( i + 0.5 ) / nchi;
        double cg = cb*ca + sb*sa*std::cos(chi);
        cg = NC::ncclamp(cg,-1.0,1.0);
        const double cg2 = cg*cg;
        const double x = R*R*cg2 + (1.0-cg2)/R;
        sum += 1.0 / ( x * std::sqrt(x) );    // = x^-1.5
      }
      return sum / nchi;
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
             || ! (preferred_orientation1.mag()>0) || ! (R1>=1e-3 && R1<=1e3) || ! (f1>0.0) )
          NCRYSTAL_THROW2( BadInput,"Invalid values in first Texture line (R must be in [1e-3,1e3])." );
      }
      else if ( texture_lines_found == 2 ) {
        if ( ! NC::safe_str2dbl( line.at(1), preferred_orientation2.at(0) )
             || ! NC::safe_str2dbl( line.at(2), preferred_orientation2.at(1) )
             || ! NC::safe_str2dbl( line.at(3), preferred_orientation2.at(2) )
             || ! NC::safe_str2dbl( line.at(4), R2 )
             || ! NC::safe_str2dbl( line.at(5), f2 )
             || ! (preferred_orientation2.mag()>0) || ! (R2>=1e-3 && R2<=1e3) || ! (f2>0.0) )
          NCRYSTAL_THROW2( BadInput,"Invalid values in second Texture line (R must be in [1e-3,1e3])." );
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

  //Precompute the per-plane cos(axis,normal) (#4) and tabulate the f-weighted texture
  //factor P(sin_theta) on a uniform grid (#1), so calcCrossSection / sampling are table
  //lookups instead of per-call azimuthal integrals.
  if ( m_has_texture ) {
    for ( auto& e : m_hklPlanes ) {
      e.cosA1 = cosAngleVectors( m_preferred_orientation1, e.hkl );
      e.cosA2 = cosAngleVectors( m_preferred_orientation2, e.hkl );
      e.texP.resize( NSINT );
      for ( int j=0; j<NSINT; ++j ) {
        const double st = double(j) / double( NSINT - 1 );
        e.texP[j] = m_f1 * sato_kernel( e.cosA1, st, m_R1 )
                  + m_f2 * sato_kernel( e.cosA2, st, m_R2 );
      }
    }
    //2-D pole-density tables for the oriented path (#6): md_pole_density(cosBetaQ,cosAlpha)
    //over cosBetaQ in [-1,1] (NBQ) and cosAlpha in [0,1] (NAL).
    m_poleTab1.resize( NBQ * NAL );
    if ( m_f2 > 0.0 ) m_poleTab2.resize( NBQ * NAL );
    for ( int i=0; i<NBQ; ++i ) {
      const double cbq = -1.0 + 2.0 * i / double( NBQ - 1 );
      for ( int j=0; j<NAL; ++j ) {
        const double ca = double(j) / double( NAL - 1 );
        m_poleTab1[ i*NAL + j ] = md_pole_density( cbq, ca, m_R1 );
        if ( m_f2 > 0.0 ) m_poleTab2[ i*NAL + j ] = md_pole_density( cbq, ca, m_R2 );
      }
    }
  }

  //Build the non-oriented total cross-section table (must be last: needs the per-plane
  //texP / cosA data above). Only the isotropic calcCrossSection uses it; the oriented
  //path is direction-dependent and keeps its per-plane evaluation.
  buildXSTable();
}

double NCP::CrystallineExtinction::poleTabInterp( const std::vector<double>& tab,
                                                  double cosBetaQ, double cosAlpha ) const
{
  //Bilinear interpolation over (cosBetaQ in [-1,1], cosAlpha in [0,1]).
  const double fb = NC::ncclamp( ( cosBetaQ + 1.0 ) * 0.5, 0.0, 1.0 ) * ( NBQ - 1 );
  const double fa = NC::ncclamp( cosAlpha, 0.0, 1.0 ) * ( NAL - 1 );
  int ib = static_cast<int>(fb); if ( ib >= NBQ-1 ) ib = NBQ-2;
  int ia = static_cast<int>(fa); if ( ia >= NAL-1 ) ia = NAL-2;
  const double tb = fb - ib, ta = fa - ia;
  const double v00 = tab[ ib*NAL + ia ],     v01 = tab[ ib*NAL + ia+1 ];
  const double v10 = tab[ (ib+1)*NAL + ia ], v11 = tab[ (ib+1)*NAL + ia+1 ];
  return ( v00*(1-ta) + v01*ta )*(1-tb) + ( v10*(1-ta) + v11*ta )*tb;
}

double NCP::CrystallineExtinction::textureFactorTab( const HKLPlane& e, double wl ) const
{
  //Linear interpolation of the precomputed P(sin_theta) table.
  const double st = 0.5 * wl / e.d_hkl;
  if ( st >= 1.0 ) return e.texP.empty() ? 1.0 : e.texP.back();
  if ( st <= 0.0 ) return e.texP.empty() ? 1.0 : e.texP.front();
  const double f = st * ( NSINT - 1 );
  const int j = static_cast<int>(f);
  const double t = f - j;
  return e.texP[j] * ( 1.0 - t ) + e.texP[j+1] * t;
}

double NCP::CrystallineExtinction::calcCrossSectionExact( double neutron_ekin ) const {

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

    const double texture_factor = m_has_texture ? textureFactorTab( e, wl ) : 1.0;
    xs_in_barns += e.strength * E_hkl * texture_factor;
  }
  xs_in_barns *= 2. * wlsq;

  return xs_in_barns;
}

void NCP::CrystallineExtinction::buildXSTable()
{
  //Build the edge-aware g(E)=xs*E grid used by the tabulated calcCrossSection (see header).
  m_xsTabE.clear();
  m_xsTabG.clear();
  if ( m_hklPlanes.empty() )
    return;

  //Bragg edge energies (neutron_ekin at wl = 2*d_hkl). One per distinct d-spacing.
  std::vector<double> edges;
  edges.reserve( m_hklPlanes.size() );
  for ( const auto& e : m_hklPlanes )
    edges.push_back( NC::wl2ekin( 2.0 * e.d_hkl ) );
  std::sort( edges.begin(), edges.end() );
  edges.erase( std::unique( edges.begin(), edges.end(),
                            []( double a, double b ){ return std::abs(a-b) <= 1e-9*std::abs(b); } ),
               edges.end() );

  const double Elo = edges.front();
  //Cover the whole epithermal tail; beyond Ehi the cross section is in its asymptotic 1/E
  //regime (all planes active, texture factor -> forward value) and is extrapolated as g/E.
  const double Ehi = std::max( edges.back() * 1000.0, 1.0e4 );

  //Base log-spaced grid (dense enough that linear interpolation of the smooth g between
  //Bragg edges is essentially exact) merged with two straddle points around every edge so
  //the Bragg step is represented sharply (lower point excludes the opening plane, upper
  //includes it).
  std::vector<double> grid;
  const double llo = std::log( Elo * ( 1.0 - 1e-6 ) );
  const double lhi = std::log( Ehi );
  const int ppd   = 150;                                   //grid points per decade
  const int nbase = std::max( 64, int( ppd * ( lhi - llo ) / std::log(10.0) ) );
  grid.reserve( nbase + 2 * edges.size() + 2 );
  for ( int k = 0; k <= nbase; ++k )
    grid.push_back( std::exp( llo + ( lhi - llo ) * k / nbase ) );
  //Around every Bragg edge add a point just below (the opening plane still closed) and a
  //geometric ladder just above it. The newly-opened reflection is at near-backscatter,
  //where extinction is strongest, so with extinction the cross section rises steeply over
  //the first few percent above the edge and needs fine post-edge sampling for linear
  //interpolation to stay accurate (without extinction this region is smooth and the extra
  //points are simply harmless).
  for ( double Ee : edges ) {
    grid.push_back( Ee * ( 1.0 - 1e-7 ) );
    grid.push_back( Ee * ( 1.0 + 1e-7 ) );
    for ( double delta = 1e-5; delta < 0.08; delta *= 1.5 )
      grid.push_back( Ee * ( 1.0 + delta ) );
  }
  std::sort( grid.begin(), grid.end() );
  grid.erase( std::unique( grid.begin(), grid.end(),
                           []( double a, double b ){ return std::abs(a-b) <= 1e-12*std::abs(b); } ),
              grid.end() );

  m_xsTabE.reserve( grid.size() );
  m_xsTabG.reserve( grid.size() );
  for ( double E : grid ) {
    m_xsTabE.push_back( E );
    m_xsTabG.push_back( calcCrossSectionExact( E ) * E );
  }

  //calcCrossSection extrapolates above the grid as g_back/ekin, which assumes g has reached
  //its asymptotic plateau at the top of the grid (above the highest Bragg edge all planes are
  //active and the texture/extinction factors sit at their near-forward values, so g=xs*E is
  //nearly constant and xs ~ 1/E). Ehi is set far above the highest edge so this holds with
  //wide margin; assert the table top is flat to ~2% to make the invariant explicit and catch
  //pathological correction parameters (assert is a no-op in release builds).
  if ( m_xsTabG.size() >= 2 ) {
    const double g_top = m_xsTabG.back();
    const double g_prev = m_xsTabG[ m_xsTabG.size() - 2 ];
    nc_assert( g_prev > 0.0 && std::abs( g_top / g_prev - 1.0 ) < 0.02 );
  }
}

double NCP::CrystallineExtinction::calcCrossSection( double neutron_ekin ) const {
  //Tabulated isotropic cross section: binary search + linear interpolation of g(E)=xs*E on
  //the edge-aware grid, then divide by E. Falls back to the exact per-plane sum when the
  //table is unavailable (e.g. no planes). The exact path is also used directly by
  //sampleScatteringEvent (which needs the per-plane weights).
  if ( m_xsTabE.empty() )
    return calcCrossSectionExact( neutron_ekin );
  if ( neutron_ekin <= m_xsTabE.front() )
    return 0.0;                                    //below the first Bragg edge: no coherent-elastic
  if ( neutron_ekin >= m_xsTabE.back() )
    return m_xsTabG.back() / neutron_ekin;         //above the tabulated range: asymptotic 1/E
  auto it = std::upper_bound( m_xsTabE.begin(), m_xsTabE.end(), neutron_ekin );
  const std::size_t i = static_cast<std::size_t>( ( it - m_xsTabE.begin() ) - 1 );
  const double e0 = m_xsTabE[i], e1 = m_xsTabE[i+1];
  const double t  = ( neutron_ekin - e0 ) / ( e1 - e0 );
  const double g  = m_xsTabG[i] * ( 1.0 - t ) + m_xsTabG[i+1] * t;
  return g / neutron_ekin;
}

NCP::CrystallineExtinction::ScatEvent NCP::CrystallineExtinction::sampleScatteringEvent( NC::RNG& rng, double neutron_ekin ) const {

  // NOTE: This sampling does not include the effect of texture and thus this plugin should not be used for Monte-Carlo simulations if texture is present.

  ScatEvent result;
  result.ekin_final = neutron_ekin;
  result.mu = 1.0;   //default (forward); overwritten once a plane is selected below

  const double wl   = NC::ekin2wl(neutron_ekin);
  const double wlsq = NC::ncsquare(wl);
  //Use the EXACT per-plane cross section here (not the tabulated calcCrossSection): the
  //inverse-CDF plane selection below accumulates the exact per-plane weights and must
  //normalise by their exact sum, otherwise the cumulative could over/undershoot 1.
  const double xs   = calcCrossSectionExact( neutron_ekin ) / ( 2. * wlsq );
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

    const double texture_factor = m_has_texture ? textureFactorTab( e, wl ) : 1.0;
    right_bound += e.strength * E_hkl * texture_factor / xs;
    nc_assert( left_bound < right_bound && right_bound <= 1.0 + 1e-9 );

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

// ************************************************************************** //
//   Anisotropic (oriented) textured coherent-elastic path (model B)         //
// ************************************************************************** //

namespace {
  //Orthonormal basis (e1,e2) perpendicular to unit vector d:
  void perpBasis( const NC::Vector& d, NC::Vector& e1, NC::Vector& e2 )
  {
    NC::Vector a = ( NC::ncabs(d[0]) < 0.9 ) ? NC::Vector(1.,0.,0.) : NC::Vector(0.,1.,0.);
    e1 = d.cross(a); e1 *= 1.0 / e1.mag();
    e2 = d.cross(e1);
  }
}

void NCP::CrystallineExtinction::setOrientation( const NCrystal::RotMatrix& crystal2lab )
{
  //Texture axes (given in the crystal/ortho frame, as the hkl normals are) -> lab:
  if ( m_has_texture ) {
    m_axis1_lab = crystal2lab * m_preferred_orientation1; m_axis1_lab *= 1.0/m_axis1_lab.mag();
    m_axis2_lab = crystal2lab * m_preferred_orientation2; m_axis2_lab *= 1.0/m_axis2_lab.mag();
  }
  m_normal_lab.clear();
  m_normal_lab.reserve( m_hklPlanes.size() );
  for ( auto& e : m_hklPlanes ) {
    NC::Vector n = crystal2lab * e.hkl;   // e.hkl is the unit reflection normal (ortho crystal frame)
    n *= 1.0 / n.mag();
    m_normal_lab.push_back( n );
  }
  m_oriented = true;
}

double NCP::CrystallineExtinction::poleDensityAtQ( std::size_t i, const NC::Vector& Qhat ) const
{
  //f-weighted modified-March-Dollase pole density of plane i at scattering-vector
  //direction Qhat (both axes & normal in lab frame):
  const NC::Vector& n = m_normal_lab[i];
  double P = m_f1 * poleTabInterp( m_poleTab1, Qhat.dot(m_axis1_lab), NC::ncabs(n.dot(m_axis1_lab)) );
  if ( m_f2 > 0.0 )
    P += m_f2 * poleTabInterp( m_poleTab2, Qhat.dot(m_axis2_lab), NC::ncabs(n.dot(m_axis2_lab)) );
  return P;
}

double NCP::CrystallineExtinction::textureFactorDir( std::size_t i, const NC::Vector& indir,
                                                     double wl ) const
{
  //Cone-average of the pole density over the Debye cone of Qhat about indir:
  const double d = m_hklPlanes[i].d_hkl;
  const double sinT = 0.5 * wl / d;
  if ( sinT > 1.0 ) return 0.0;
  const double cosT = std::sqrt( std::max(0.0, 1.0 - sinT*sinT) );
  NC::Vector e1, e2; perpBasis( indir, e1, e2 );
  const unsigned npsi = 64;
  double sum = 0.0;
  for ( unsigned j=0; j<npsi; ++j ) {
    const double psi = NC::k2Pi * ( j + 0.5 ) / npsi;
    //Qhat(psi) = sinT*indir + cosT*(cos psi e1 + sin psi e2)   (Qhat.indir = sinT)
    NC::Vector Q = indir * sinT + e1 * (cosT*std::cos(psi)) + e2 * (cosT*std::sin(psi));
    sum += poleDensityAtQ( i, Q );
  }
  return sum / npsi;
}

double NCP::CrystallineExtinction::calcCrossSectionDir( double neutron_ekin,
                                                        const NC::Vector& indir ) const
{
  if ( !m_oriented || !m_has_texture )
    return calcCrossSection( neutron_ekin );   //no texture frame -> isotropic value

  double xs_in_barns = 0.0;
  const double wl   = NC::ekin2wl( neutron_ekin );
  const double wlsq = NC::ncsquare( wl );
  const double mu0  = 0.;
  for ( std::size_t i=0; i<m_hklPlanes.size(); ++i ) {
    auto& e = m_hklPlanes[i];
    if ( wl > 2 * e.d_hkl ) break;
    double E_hkl = 1.0;
    if ( m_has_extinction ) {
      if      ( m_model_option == 0 ) E_hkl = uncorr_blk_mdl( m_Nc, wl, e.F_hkl, m_l, e.d_hkl, mu0, m_Gg, m_L, m_tilt_dist_option );
      else if ( m_model_option == 1 ) E_hkl = corr_blk_mdl( m_Nc, wl, e.F_hkl, m_l, e.d_hkl, mu0, m_Gg, m_L );
      else if ( m_model_option == 2 ) E_hkl = BC_pure_extn_mdl( m_Nc, wl, e.F_hkl, m_l, e.d_hkl, m_Gg, m_L, m_tilt_dist_option, m_recipe );
      else if ( m_model_option == 3 ) E_hkl = BC_mix_extn_mdl( m_Nc, wl, e.F_hkl, m_l, e.d_hkl, m_Gg, m_L, m_tilt_dist_option, m_recipe );
      else                            E_hkl = BC_mod_extn_mdl( m_Nc, wl, e.F_hkl, m_l, e.d_hkl, m_Gg, m_L, m_tilt_dist_option, m_recipe );
    }
    xs_in_barns += e.strength * E_hkl * textureFactorDir( i, indir, wl );
  }
  return xs_in_barns * 2.0 * wlsq;
}

NCP::CrystallineExtinction::ScatEventDir
NCP::CrystallineExtinction::sampleScatteringEventDir( NC::RNG& rng, double neutron_ekin,
                                                      const NC::Vector& indir ) const
{
  ScatEventDir result;
  result.ekin_final = neutron_ekin;   //elastic
  result.outdir = indir;

  if ( !m_oriented || !m_has_texture ) {
    //Fall back to isotropic mu sampling, expanded to a random azimuth:
    auto ev = sampleScatteringEvent( rng, neutron_ekin );
    NC::Vector e1, e2; perpBasis( indir, e1, e2 );
    const double phi = NC::k2Pi * rng.generate();
    const double smu = std::sqrt( std::max(0.0, 1.0 - ev.mu*ev.mu) );
    result.outdir = indir * ev.mu + e1 * (smu*std::cos(phi)) + e2 * (smu*std::sin(phi));
    return result;
  }

  const double wl = NC::ekin2wl( neutron_ekin );
  const double mu0 = 0.;

  //1) choose plane i with probability proportional to its directional partial xs:
  std::vector<double> cw; cw.reserve( m_hklPlanes.size() );
  std::vector<std::size_t> idx; idx.reserve( m_hklPlanes.size() );
  double acc = 0.0;
  for ( std::size_t i=0; i<m_hklPlanes.size(); ++i ) {
    auto& e = m_hklPlanes[i];
    if ( wl > 2 * e.d_hkl ) break;
    double E_hkl = 1.0;
    if ( m_has_extinction ) {
      if      ( m_model_option == 0 ) E_hkl = uncorr_blk_mdl( m_Nc, wl, e.F_hkl, m_l, e.d_hkl, mu0, m_Gg, m_L, m_tilt_dist_option );
      else if ( m_model_option == 1 ) E_hkl = corr_blk_mdl( m_Nc, wl, e.F_hkl, m_l, e.d_hkl, mu0, m_Gg, m_L );
      else if ( m_model_option == 2 ) E_hkl = BC_pure_extn_mdl( m_Nc, wl, e.F_hkl, m_l, e.d_hkl, m_Gg, m_L, m_tilt_dist_option, m_recipe );
      else if ( m_model_option == 3 ) E_hkl = BC_mix_extn_mdl( m_Nc, wl, e.F_hkl, m_l, e.d_hkl, m_Gg, m_L, m_tilt_dist_option, m_recipe );
      else                            E_hkl = BC_mod_extn_mdl( m_Nc, wl, e.F_hkl, m_l, e.d_hkl, m_Gg, m_L, m_tilt_dist_option, m_recipe );
    }
    acc += e.strength * E_hkl * textureFactorDir( i, indir, wl );
    cw.push_back( acc ); idx.push_back( i );
  }
  if ( acc <= 0.0 ) return result;   //no scattering (shouldn't happen if xs>0)

  const double r = rng.generate() * acc;
  std::size_t pick = idx.back();
  for ( std::size_t j=0; j<cw.size(); ++j ) { if ( r < cw[j] ) { pick = idx[j]; break; } }

  //2) sample azimuth psi on the Debye cone with density ~ poleDensityAtQ(pick,Qhat(psi)).
  //   Inverse-CDF (histogram) sampler: guaranteed to terminate and uses no rejected
  //   sample, robust for any (finite) texture strength.
  const double d = m_hklPlanes[pick].d_hkl;
  const double sinT = 0.5 * wl / d;
  const double cosT = std::sqrt( std::max(0.0, 1.0 - sinT*sinT) );
  NC::Vector e1, e2; perpBasis( indir, e1, e2 );
  auto Qof = [&]( double psi ) {
    return indir * sinT + e1 * (cosT*std::cos(psi)) + e2 * (cosT*std::sin(psi));
  };
  constexpr int NPSI = 128;
  double cdf[NPSI]; double s2 = 0.0;
  for ( int j=0; j<NPSI; ++j ) {
    s2 += poleDensityAtQ( pick, Qof( NC::k2Pi * ( j + 0.5 ) / NPSI ) );
    cdf[j] = s2;
  }
  double psi;
  if ( s2 <= 0.0 ) {
    psi = NC::k2Pi * rng.generate();   //degenerate (should not happen) -> uniform
  } else {
    const double tgt = rng.generate() * s2;
    int jb = 0; while ( jb < NPSI-1 && tgt > cdf[jb] ) ++jb;
    psi = NC::k2Pi * ( jb + rng.generate() ) / NPSI;   //uniform within the chosen bin
  }
  NC::Vector Q = Qof( psi );
  //3) elastic outgoing direction: k_hat' = k_hat - 2 sinT * Qhat
  result.outdir = indir - Q * ( 2.0 * sinT );
  result.outdir *= 1.0 / result.outdir.mag();
  return result;
}
