//
//  TPZHdivTransfer.h
//
//  Created by Jorge Paúl Ordóñez Andrade on 15/07/19.
//

#ifndef TPZHdivTransfer_h
#define TPZHdivTransfer_h

#include <stdio.h>
#include "pzvec.h"
#include "pzmatrix.h"

template<class TVar>
class TPZHdivTransfer : public TPZMatrix<TVar> {
    
    TPZVec<int64_t> fIndexes;
    
public:
    
    /** @brief Default constructor */
    TPZHdivTransfer();
    
    /** @brief Copy constructor */
    TPZHdivTransfer(const TPZHdivTransfer<TVar> &cp);
    
    CLONEDEF(TPZHdivTransfer<TVar>)
    TPZHdivTransfer* NewMatrix() const override{ return new TPZHdivTransfer{};}

    /** @brief Creates a copy from another TPZBlockDiagonal*/
    void CopyFrom(const TPZMatrix<TVar> *  mat) override
    {                                                           
        auto *from = dynamic_cast<const TPZHdivTransfer<TVar> *>(mat);                
        if (from) {                                               
        *this = *from;                                          
        }                                                         
        else                                                      
        {                                                       
            PZError<<__PRETTY_FUNCTION__;                         
            PZError<<"\nERROR: Called with incompatible type\n."; 
            PZError<<"Aborting...\n";                             
            DebugStop();                                          
        }                                                       
    }
    
    /** @brief Constructor based on indexes */
    TPZHdivTransfer(int64_t rows, int64_t cols, TPZVec<int64_t> &Indexes);
    
    /** @brief Default constructor */
    ~TPZHdivTransfer();
    
    /**
     * @brief It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
     * @param x Is x on the above operation
     * @param y Is y on the above operation
     * @param z Is z on the above operation
     * @param alpha Is alpha on the above operation
     * @param beta Is beta on the above operation
     * @param opt Indicates if is Transpose or not
     */
    virtual void MultAdd(const TPZFMatrix<TVar> & x,const TPZFMatrix<TVar>& y, TPZFMatrix<TVar>& z,
                         const TVar alpha=1., const TVar beta = 0., const int opt = 0) const override;
    

    /** @brief Copy constructor */
    virtual void Print(const char *name = NULL, std::ostream &out = std::cout , const MatrixOutputFormat form = EFormatted) const override;
    
    /** @brief Gather y elements into x vector */
    void Gather(TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &x);
    
    /** @brief Scatter y elements into x vector */
    void Scatter(TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &x);
    
    /** @brief Set indexes */
    void SetIndexes(int64_t rows, int64_t cols, TPZVec<int64_t> &Indexes);
    
    
    /** @brief Get indexes */
    TPZVec<int64_t> & GetIndexes();

protected:
	inline TVar *&Elem() override
    {
        DebugStop();//Put the values of TPZMatrix or call the TPZMatrix function
        // return fIndexes.begin();
    }
    inline const TVar *Elem() const override
    {
        DebugStop();
        // return fIndexes.begin();
    }
    inline int64_t Size() const override
    {
        DebugStop();
        return fIndexes.size();
    }
};

#endif /* TPZHdivTransfer_h */
