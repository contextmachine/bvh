//
// Created by Andrew Astakhov on 23.01.25.
//

#ifndef ACOUSTICS_H
#define ACOUSTICS_H
#include "vec.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <math.h>
#define DEBUG 1
namespace acoustics {
    template<class T>
    // Сommon material absorption coefficients $ alpha_0 $
    struct  MaterialAbsorption{
        static constexpr T concrete=0.03;
        static constexpr T carpetedFloor=0.4;
        static constexpr T curtains=0.5; // шторы
        static constexpr T audience=0.6; // люди блять
    };
    template<class T>
    // Сommon absorption coefficients $ alpha_0 $ for a interior elements
    struct  Absorption{
        static constexpr T walls=0.02;
        static constexpr T floor=0.4;
        static constexpr T ceiling=0.3;

    };


    template<class T>
    // P: Source sound power in watts.
    //Reference distance  (commonly 1 meter). НЕ ПУТАТЬ С РЕАЛЬНЫМ РАССТОЯНИЕМ КОТОРОЕ МЫ СОБИРАЕМСЯ ВЫЧИСЛЯТЬ
    // returns  $ I_0 $ $ W/m^2 $
    constexpr T initialIntensity(const T P, T r0=1) {
        return (P/(4*M_PI*std::pow(r0,2)));
    }


    template<class T>
    // I_{\text{ref}} = 10^{-12} \, \text{W/m}^2  (reference intensity for air).
    static constexpr T I_ref= std::pow<T>(10.,-12);

    template<class T>
    // Initial Sound Pressure Level coefficients $ SPL_0 $
    constexpr T initialSoundPressure(const T P, const T r0=1) {
        return 10*log10(initialIntensity(P,r0)/I_ref<T>);
    }
    template<class T>
    struct SoundSource {
        T power;
        T r0=1;
        T SPL_0;
        SoundSource()=default;
        SoundSource(T P, T r0):power(P),r0(r0),SPL_0(initialSoundPressure(power,r0)){}
        SoundSource(T P):power(P),r0(1.),SPL_0(initialSoundPressure(power,r0)){}

    };

    template<class T>
    /*$ alpha_{air} $ Аir absorption coefficient (in dB/m), frequency-dependent.
    Typical values for  $ alpha_{air} $  (500–1000 Hz, 20°C, 50% RH): */
    struct  AirAbsorption{
        static constexpr T lowFreq=0.01;// dB/m for low frequencies.
        static constexpr T highFreq=0.03;// dB/m for high frequencies.

    };

    template<class T>
    T reflectionLoss(const bvh::vec<T,3> &d, const bvh::vec<T,3> &n,const T alpha_0){
        #ifdef DEBUG
        assert( bvh::detail::is_close<T>(d.length(),1.), "reflectionLoss: ray direction should be normalized!");
        assert( bvh::detail::is_close<T>(n.length(),1.), "reflectionLoss: normal should be normalized!");
        #endif
        return  alpha_0*d.dot(n);
    }
    template<class T>
    T airAbsorptionLoss(const T distance,const T alpha_air=AirAbsorption<T>::highFreq) {
        return std::exp(-alpha_air * distance);

    }

    template<class T>
    constexpr T distanceLoss(const T distance,const SoundSource<T> &sound_source) {

        if ( bvh::detail::is_close<T>(distance,0.)) {
            std::cerr<<"distanceLoss: distance="<<distance<<"\n";

            return std::numeric_limits<T>::infinity();
        }
        return sound_source.SPL_0 - 20 * std::log10(distance);

    }

    template<class T>
    constexpr T combineLoss(
        const T distance,
        const T reflection_loss,
        const T initial_intensity,
        const T alpha_air=AirAbsorption<T>::highFreq) {
        auto air_absorption_loss=airAbsorptionLoss(distance,alpha_air);
        return initial_intensity * (1 - reflection_loss) * air_absorption_loss / std::pow(distance , 2);
    }
    template<class T>
    /*
     *
\phi(t) = \pi t \quad (0 \leq t \leq 1)

     */
    class Loxodrome {
    public:
        T phi(const T t) {
            return t*M_PI;
        }




    };

    template<class T>
    class AcousticsSimulator {
        public:

        AcousticsSimulator()=default;
        ~AcousticsSimulator()=default;


    };



}

#endif //ACOUSTICS_H
