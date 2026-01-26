//
// Created by Maria on 22.03.2024 kolupanova@inr.ru
//

#ifndef UNPACKING_HIT_H
#define UNPACKING_HIT_H

#include <iostream>
#include <cmath>
#include <vector>
#include "TObject.h"

struct ChannelData : public TObject{
    bool _channels_flag{true};
    unsigned int _channel_id;
    double _time{NAN};
    double _quality{NAN};

    double _amplitude{NAN};
    double _chargeADC{NAN};
    double _chargePE_Integral{NAN};
    double _chargePE_Amplitude{NAN};
    double _chargePE_Integral_XTalk_with_background{NAN};
    double _chargePE_Amplitude_XTalk_with_background{NAN};
    double _chargePE_Integral_XTalk_without_background{NAN};
    double _chargePE_Amplitude_XTalk_without_background{NAN};

    double _reconstructed_amplitude{NAN};
    double _reconstructed_chargePE_Amplitude{NAN};
    int _overflow;

    ClassDef(ChannelData, 1);
};

struct PrototypeHit : public TObject  {
    double _x{NAN};
    double _y{NAN};
    double _Ax{NAN};
    double _Bx{NAN};
    double _Ay{NAN};
    double _By{NAN};
    int _quality;
    std::vector<ChannelData> _channelsData;

    ClassDef(PrototypeHit, 1);
};

#endif //UNPACKING_HIT_H
