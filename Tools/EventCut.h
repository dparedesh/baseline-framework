#ifndef EVENTCUT_H
#define EVENTCUT_H


#include "TString.h"


class EventCut {

public:

    EventCut(TString cut,TString name);
    ~EventCut();
 
    inline TString GetName(){return m_name;};
    inline TString GetCut(){return m_cut;};
    inline void SetMinimum(float min){m_min=min;};
    inline void SetMaximum(float max){m_max=max;};
    inline float GetMinimum(){return m_min;};
    inline float GetMaximum(){return m_max;};
    inline void SetCut(TString cut){m_cut=cut;};
    


private:

    TString m_cut;
    TString m_name;
    float m_min;
    float m_max;

};

#endif // EVENTCUT_J
