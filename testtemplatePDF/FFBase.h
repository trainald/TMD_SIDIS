//
//  FFBase.h
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 10/24/22.
//

#ifndef FFBase_h
#define FFBase_h

class FFBase
{
public:
    
    // Constructors -----------------------------------------------------------
    
    FFBase()
    {};
    
    FFBase(double Q, double muQ, int hadron, int flavour)// changed order hadron flavor
    :   m_HardScale(Q),
        m_RenormalizationScale(muQ),
        m_Flavour(flavour),
        m_Hadron(hadron)
    {};
    
    // Methods - Get ----------------------------------------------------------
        
    double Get_HardScale() const;
    
    double Get_RenormalizationScale() const;
    
    int Get_Flavour() const;
    
    int Get_Hadron() const;
    
    
    // Methods - Set ----------------------------------------------------------
    
    void Set_HardScale(double Q);
    
    void Set_RenormalizationScale(double muQ);
    
    void Set_Flavour(int flavour);
    
    void Set_Hadron(int hadron);

protected:
    
    // Member variables -------------------------------------------------------
    
    /* \brief hard scale.*/
    double m_HardScale;
    
    /* \brief renormalization scale.*/
    double m_RenormalizationScale;
    
    /* \brief parton's flavour.*/
    int m_Flavour;
    
    /* \brief probed hadron.*/
    int m_Hadron;
};


// Definition -------------------------------------------------------------

// Get -----------------------------------------------------
// Get Hard Scale Q
double FFBase::Get_HardScale() const
{
    return m_HardScale;
}

// Get Renormalization Scale muQ
double FFBase::Get_RenormalizationScale() const
{
    return m_RenormalizationScale;
}

// Get parton's flavour
int FFBase::Get_Flavour() const
{
    return m_Flavour;
}

// Get probed hadron
int FFBase::Get_Hadron() const
{
    return m_Hadron;
}

// Set -----------------------------------------------------
// Set Hard Scale Q
void FFBase::Set_HardScale(double Q)
{
    m_HardScale = Q;
}

// Set Renormalization Scale muQ
void FFBase::Set_RenormalizationScale(double muQ)
{
    m_RenormalizationScale = muQ;
}

// Set parton's flavour
void FFBase::Set_Flavour(int flavour)
{
    m_Flavour = flavour;
}

// Set probed hadron
void FFBase::Set_Hadron(int hadron)
{
    m_Hadron = hadron;
}

#endif /* FFBase_h */
