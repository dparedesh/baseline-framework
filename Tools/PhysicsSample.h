#ifndef PHYSICSSAMPLE_H
#define PHYSICSSAMPLE_H


#include "TString.h"


class PhysicsSample {

public:

    PhysicsSample(TString pFile, double xs);
    ~PhysicsSample();
 
    inline TString GetFile(){return m_pfile;};
    inline double GetXS(){return m_xs;};
    inline double GetScaleFactor(){return m_sf;};

    inline void SetScaleFactor(float sf){m_sf=sf;};


private:

    TString m_pfile;
    double m_xs;
    double m_sf;

};

#endif // PHYSICSSAMPLE_H
