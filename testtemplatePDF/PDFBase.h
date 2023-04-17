//
//  PDFBase.h
//  testtemplatePDF
//
//  Created by Tommaso Rainaldi on 10/21/22.
//

#ifndef PDFBase_h
#define PDFBase_h

class PDFBase
{
public:
    
    // Constructors -----------------------------------------------------------
    
    PDFBase()
    {};
    
    PDFBase(double Q, double muQ, int flavour, int hadron)
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
double PDFBase::Get_HardScale() const
{
    return m_HardScale;
}

// Get Renormalization Scale muQ
double PDFBase::Get_RenormalizationScale() const
{
    return m_RenormalizationScale;
}

// Get parton's flavour
int PDFBase::Get_Flavour() const
{
    return m_Flavour;
}

// Get probed hadron
int PDFBase::Get_Hadron() const
{
    return m_Hadron;
}

// Set -----------------------------------------------------
// Set Hard Scale Q
void PDFBase::Set_HardScale(double Q)
{
    m_HardScale = Q;
}

// Set Renormalization Scale muQ
void PDFBase::Set_RenormalizationScale(double muQ)
{
    m_RenormalizationScale = muQ;
}

// Set parton's flavour
void PDFBase::Set_Flavour(int flavour)
{
    m_Flavour = flavour;
}

// Set probed hadron
void PDFBase::Set_Hadron(int hadron)
{
    m_Hadron = hadron;
}


#endif /* PDFBase_h */
