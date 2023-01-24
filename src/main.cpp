#include <iostream>
#include <memory>
#include <string.h>
#include <math.h>

#include "StochBlock_types.h"


template < typename T >
class Borders {

public:

    Borders(){}
    Borders( const size_t rows, const size_t cols, const size_t slices = 0 );
    Borders( const T & p_lower, const T & p_upper ) : m_lower( p_lower ), m_upper( p_upper ) {}
//    virtual ~Borders() {}

    T getLower() const { return m_lower; }
    T getUpper() const { return m_upper; }

    double getLowerAt( const size_t i, const size_t j, const size_t r = 0 ) const { return m_lower.at( i, j ); }
    double getUpperAt( const size_t i, const size_t j, const size_t r = 0 ) const { return m_upper.at( i, j ); }

protected:

    T m_lower;
    T m_upper;

};

template < typename T >
Borders<T>::Borders( const size_t rows, const size_t cols, const size_t /* slices */ )
{
    m_lower = T( rows, cols );
    m_upper = T( rows, cols );

    m_lower.fill( R_NegInf );
    m_upper.fill( R_PosInf );
}

template <>
Borders<Array>::Borders( const size_t rows, const size_t cols, const size_t slices )
{
    m_lower = Array( rows, cols, slices );
    m_upper = Array( rows, cols, slices );

    m_lower.fill( R_NegInf );
    m_upper.fill( R_PosInf );
}

template <>
double Borders<Array>::getLowerAt( const size_t i, const size_t j, const size_t r ) const
{
    return m_lower.at( i, j, r );
}

template <>
double Borders<Array>::getUpperAt( const size_t i, const size_t j, const size_t r ) const
{
    return m_upper.at( i, j, r );
}




// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

// Exposed functions
// [[Rcpp::export(name = ".meanByBlocks")]]
Rcpp::List meanByBlocks( const Array & M, const IVector & clu, const IVector & nClu, const IVector & n, const std::string diagonal = "ignore",
                         const std::string & sBorders = "none", const Rcpp::Nullable<Array> & bordersMatLower = R_NilValue, const Rcpp::Nullable<Array> & bordersMatUpper = R_NilValue,
                         const Rcpp::Nullable<DMatrix> & bordersSeperateLower = R_NilValue, const Rcpp::Nullable<DMatrix> & bordersSeperateUpper = R_NilValue, const bool & addOne = true, const double & eps = 0.001);
// [[Rcpp::export(name = ".kmBlock")]]
Rcpp::List kmBlock( const Array & M, const IVector & clu, const Array & weights, const DVector & uWeights, const IVector & n, const IVector & nClu, const std::string & diagonal = "ignore", const double weightClusterSize=1.0, 
                    const std::string & sBorders = "none", const Rcpp::Nullable<Array> & bordersMatLower = R_NilValue, const Rcpp::Nullable<Array> & bordersMatUpper = R_NilValue,
                    const Rcpp::Nullable<DMatrix> & bordersSeperateLower = R_NilValue, const Rcpp::Nullable<DMatrix> & bordersSeperateUpper = R_NilValue, const int & maxNoImp = 0, const bool & addOne = true, const double & eps = 0.001);
// [[Rcpp::export(name = ".critFunction")]]
double critFunction( const Array & M, const IVector & clu, const Array & weights, const DVector & uWeights, const int dimensions, const IVector & n, const double weightClusterSize = 1, const std::string & diagonal = "ignore",
                     const std::string & sBorders = "none", const Rcpp::Nullable<Array> & bordersMatLower = R_NilValue, const Rcpp::Nullable<Array> & bordersMatUpper = R_NilValue,
                     const Rcpp::Nullable<DMatrix> & bordersSeperateLower = R_NilValue, const Rcpp::Nullable<DMatrix> & bordersSeperateUpper = R_NilValue, const bool & addOne = true, const double & eps = 0.001);

// Functions forward declarations
void meansByBlocks( const Array & M, Array & res, const IVector & clu, const IVector & nClu, DMatrix & p_pSepare, const Array & superBlockMeans, const DMatrix & superBlockDiagMeans, const IVector & n, BorderType p_type,
                    const Borders<Array> & p_btBorders, const Borders<DMatrix> & p_btBordersSeperate, const Diagonale sDiagonal  = Diagonale::Ignore, const bool & addOne = true, const double & eps = 0.001);
double criterialFunction( const Array & M, const IVector & clu, const Array & weights, const DVector & uWeights, const Array & meansMat, const DMatrix & p_mSeparate, const Diagonale p_diagonale, const DVector & logProbGroups, const double weightClusterSize);

void setGroups( const Array & M, IVector & clu, const Array & weights, const DVector & uWeights, const Array & meansMat, const IVector & nClu, const IVector & n, const DMatrix & p_mSeparate, const Diagonale p_diagonale, const double weightClusterSize, const IVector & cluMode, IVector & countGroups, DVector & logProbGroups);

unsigned int belongsTo( const int & group, const IVector & borders );

DMatrix relationsMeans( const Array & M, const IVector & n );

double meanMatrix( const DMatrix & p_matrix );

void superblockMeansFun( const Array & M, const IVector & n, const Diagonale p_diagonale, const IVector & unitMode, Array & superBlockMeans, DMatrix & superBlockDiagMeans);

Diagonale getDiagonale( const std::string & p_sDiagonal );
BorderType getBorderType( const std::string & p_sBorder );

void checkInputBorders( const Diagonale & p_diagonale,
                    const Rcpp::Nullable<Array> & bordersMatLower, const Rcpp::Nullable<Array> & bordersMatUpper,
                    const Rcpp::Nullable<DMatrix> & bordersSeperateLower, const Rcpp::Nullable<DMatrix> & bordersSeperateUpper );

std::ostream & operator << ( std::ostream & p_stream, const Diagonale p_diag );
					
Rcpp::List kmBlock( const Array & M, const IVector & clu, const Array & weights, const DVector & uWeights, const IVector & n, const IVector & nClu, const std::string & diagonal, const double weightClusterSize, 
                    const std::string & sBorders, const Rcpp::Nullable<Array> & bordersMatLower, const Rcpp::Nullable<Array> & bordersMatUpper,
                    const Rcpp::Nullable<DMatrix> & bordersSeperateLower, const Rcpp::Nullable<DMatrix> & bordersSeperateUpper, const int & maxNoImp, const bool & addOne, const double & eps)
{
    const Diagonale dDiag = getDiagonale( diagonal );
    const BorderType eBorders = getBorderType( sBorders );
    const int K = Rcpp::sum( nClu );
	IVector cumN = Rcpp::cumsum( n );	
	// Rcpp::Rcout << "OK0\n";
	IVector unitMode(clu.size());

	// Rcpp::Rcout << "clu: " << clu << "\n";
	// Rcpp::Rcout << "clu.size(): " << clu.size() << "\n";
	// Rcpp::Rcout << "unitMode.size(): " << unitMode.size() << "\n";
	// Rcpp::Rcout << "unitMode: " << unitMode << "\n";
		
	int ind=0;
	for(int i = 0; i<clu.size();++i){
		if(i==cumN.at(ind)){
			++ind;
		}
		// Rcpp::Rcout << "i: " << i << ", ind: " << ind << "\n";
		unitMode.at(i)=ind;
	}
	// Rcpp::Rcout << "OK1\n";
	
	IVector countGroups( K );
	IVector cluMode( K );
	ind = 0;
	for(int i = 0;i<nClu.size();++i){
		for(int j = 0; j<nClu.at(i);++j){
			cluMode.at(ind)=i;
//			Rcpp::Rcout << "cluMode(" << ind <<") = "<< i << "\n";
			++ind;
		}
	}
	// Rcpp::Rcout << "OK2\n";
	
	for( unsigned int i = 0; i < static_cast<unsigned int>( clu.size() ); ++i ) {
		++countGroups.at(clu.at(i));
	}
	// Rcpp::Rcout << "OK3\n";	
	DVector logProbGroups( K );
	for(int i = 0; i< K;++i){
		logProbGroups.at(i)=log(static_cast<double>(countGroups.at(i))/n.at(cluMode.at(i)));
	}
	
	// Rcpp::Rcout << "OK4\n";	
    Borders<Array> bordersMeanstMat;
    Borders<DMatrix> bordersSeperate;

//    if( useBorders ) {
//        checkInputBorders( dDiag, bordersMatLower, bordersMatUpper, bordersSeperateLower, bordersSeperateUpper );
//        if( dDiag == Diagonale::Seperate ) {
//            bordersSeperate = Borders<DMatrix>( Rcpp::as<DMatrix>( bordersSeperateLower ), Rcpp::as<DMatrix>( bordersSeperateUpper ) );
//        }
//        bordersMeanstMat = Borders<Array>( Rcpp::as<Array>( bordersMatLower ), Rcpp::as<Array>( bordersMatUpper ) );
//    }
//    else {
//        if( dDiag == Diagonale::Seperate ) {
//            bordersSeperate = Borders<DMatrix>( K, M.n_slices );
//        }
//        bordersMeanstMat = Borders<Array>( K, K, M.n_slices );
//    }

    if( eBorders != BorderType::None ) {
        checkInputBorders( dDiag, bordersMatLower, bordersMatUpper, bordersSeperateLower, bordersSeperateUpper );
        if( dDiag == Diagonale::Seperate ) {
            bordersSeperate = Borders<DMatrix>( Rcpp::as<DMatrix>( bordersSeperateLower ), Rcpp::as<DMatrix>( bordersSeperateUpper ) );
        }
        bordersMeanstMat = Borders<Array>( Rcpp::as<Array>( bordersMatLower ), Rcpp::as<Array>( bordersMatUpper ) );
    }
    else {
        if( dDiag == Diagonale::Seperate ) {
            bordersSeperate = Borders<DMatrix>( K, M.n_slices );
        }
        bordersMeanstMat = Borders<Array>( K, K, M.n_slices );
    }

	// Rcpp::Rcout << "OK5\n";	
	
    DMatrix pSeparate;
    Array meanBlocks, mSeparate;


	const size_t S = static_cast<size_t>( n.size() );
	Array superBlockMeans(S, S, M.n_slices, arma::fill::zeros );
	DMatrix superBlockDiagMeans( S, M.n_slices, arma::fill::zeros );
	
	superblockMeansFun(M, n, dDiag, unitMode, superBlockMeans, superBlockDiagMeans);

	// Rcpp::Rcout << "superBlockMeans \n";
	// for(unsigned int r = 0; r<M.n_slices; ++r){
	// Rcpp::Rcout << "Slice "<< r << "\n";
		// for(unsigned int i = 0;i<S;++i){
			// for(unsigned int j = 0; j<S;++j){
				// Rcpp::Rcout << superBlockMeans(i,j,r) << "\t";
			// }
			// Rcpp::Rcout << "\n";
		// }
	// }
	// Rcpp::Rcout << "\n";


	// Rcpp::Rcout << "superBlockDiagMeans \n";
	// for(unsigned int i = 0;i<S;++i){
		// for(unsigned int r = 0; r<M.n_slices; ++r){
			// Rcpp::Rcout << superBlockDiagMeans(i,r) << "\t";
		// }
		// Rcpp::Rcout << "\n";
	// }
	// Rcpp::Rcout << "\n";
	
//    const DMatrix MEANS = relationsMeans( M, n );
	// Rcpp::Rcout << "OK6\n";	
    meansByBlocks( M, meanBlocks, clu, nClu, pSeparate, superBlockMeans, superBlockDiagMeans, n, eBorders, bordersMeanstMat, bordersSeperate, dDiag, addOne, eps);
	// Rcpp::Rcout << "OK7\n";	
	
	// Rcpp::Rcout << "cluMode \n";
	// Rcpp::Rcout << cluMode <<"\n";
	// Rcpp::Rcout << "countGroups \n";
	// Rcpp::Rcout << countGroups <<"\n";
	// Rcpp::Rcout << "logProbGroups \n";
	// Rcpp::Rcout << logProbGroups <<"\n";
	
	
	
	// Rcpp::Rcout << "superBlockMeans \n";
	// Rcpp::Rcout << superBlockMeans <<"\n";
	// Rcpp::Rcout << "superBlockDiagMeans \n";
	// Rcpp::Rcout << superBlockDiagMeans <<"\n";
	
	// Rcpp::Rcout << "clu \n";
	// Rcpp::Rcout << clu <<"\n";

	// Rcpp::Rcout << "meanBlocks \n";
	// Rcpp::Rcout << meanBlocks <<"\n";
	// Rcpp::Rcout << "\n";	

	
	
//	for(unsigned int r = 0; r<M.n_slices; ++r){
//	Rcpp::Rcout << "Slice "<< r << "\n";
//		for(unsigned int i = 0;i<K;++i){
//			for(unsigned int j = 0; j<K;++j){
//				Rcpp::Rcout << meanBlocks(i,j,r) << "\t";
//			}
//			Rcpp::Rcout << "\n";
//		}
//	}


	

    IVector newClu = Rcpp::clone( clu );
    IVector bestClu;

	// Rcpp::Rcout << "newClu \n";
	// Rcpp::Rcout << newClu <<"\n";	

    meansByBlocks( M, meanBlocks, newClu, nClu, pSeparate, superBlockMeans, superBlockDiagMeans, n, eBorders, bordersMeanstMat, bordersSeperate, dDiag, addOne, eps);
	// Rcpp::Rcout << "OK7\n";	
	// Rcpp::Rcout << "meanBlocks \n";
	// Rcpp::Rcout << meanBlocks <<"\n";
	// Rcpp::Rcout << "\n";
	
	// Rcpp::Rcout << "M\n" << M <<"\n";		
	// Rcpp::Rcout << "weights\n" << weights <<"\n";		
	// Rcpp::Rcout << "meanBlocks\n" << meanBlocks <<"\n";		
	// Rcpp::Rcout << "pSeparate\n" << pSeparate <<"\n";		
	// Rcpp::Rcout << "dDiag\n" << dDiag <<"\n";
	// Rcpp::Rcout << "logProbGroups\n" << logProbGroups <<"\n";
	// Rcpp::Rcout << "weightClusterSize\n" << weightClusterSize <<"\n";

    double newCf = criterialFunction( M, newClu, weights, uWeights, meanBlocks, pSeparate, dDiag, logProbGroups, weightClusterSize);
    double bestCf = DBL_MAX;
	// Rcpp::Rcout << "OK8\n";	

	// Rcpp::Rcout << "bestCf: " << bestCf <<"\n";
	// Rcpp::Rcout << "newCf: " << newCf <<"\n";
	bestClu = Rcpp::clone( newClu );
	bestCf = newCf;
	int noImp= 0;
    while( noImp <= maxNoImp) {
        setGroups( M, newClu, weights, uWeights, meanBlocks, nClu, n, pSeparate, dDiag, weightClusterSize, cluMode, countGroups, logProbGroups);
        meansByBlocks( M, meanBlocks, newClu, nClu, pSeparate, superBlockMeans, superBlockDiagMeans, n, eBorders, bordersMeanstMat, bordersSeperate, dDiag, addOne, eps);
		// Rcpp::Rcout << "newClu \n" << newClu <<"\n";		
		// Rcpp::Rcout << "bestClu \n" << bestClu <<"\n";		
        newCf = criterialFunction( M, newClu, weights, uWeights, meanBlocks, pSeparate, dDiag, logProbGroups, weightClusterSize);
		// Rcpp::Rcout << "bestCf: " << bestCf <<"\n";
		// Rcpp::Rcout << "newCf: " << newCf <<"\n";
		if (newCf < bestCf){
			bestClu = Rcpp::clone( newClu );
			bestCf = newCf;
			noImp = 0;
		}else{
			noImp++;
		}
    }
	
	
	meansByBlocks( M, meanBlocks, bestClu, nClu, pSeparate, superBlockMeans, superBlockDiagMeans, n, eBorders, bordersMeanstMat, bordersSeperate, dDiag, addOne, eps);

	for(int i = 0; i< K;++i){
		countGroups.at(i)=0;
	}

	for( unsigned int i = 0; i < static_cast<unsigned int>( bestClu.size() ); ++i ) {
		++countGroups.at(bestClu.at(i));
	}
	for(int i = 0; i< K;++i){
		logProbGroups.at(i)=log(static_cast<double>(countGroups.at(i))/n.at(cluMode.at(i)));
	}
	// Rcpp::Rcout << "meanBlocks\n" << meanBlocks <<"\n";		
	// Rcpp::Rcout << "pSeparate\n" << pSeparate <<"\n";		
	// Rcpp::Rcout << "dDiag\n" << dDiag <<"\n";
	// Rcpp::Rcout << "logProbGroups\n" << logProbGroups <<"\n";
	// Rcpp::Rcout << "weightClusterSize\n" << weightClusterSize <<"\n";

	// double tmpCf = criterialFunction( M, bestClu, weights, meanBlocks, pSeparate, dDiag, logProbGroups, weightClusterSize);
	// Rcpp::Rcout << "tmpCf: " << tmpCf <<"\n";
	
    return Rcpp::List::create( Rcpp::Named( "bestCf" ) = bestCf, Rcpp::Named( "bestClu" ) = bestClu, Rcpp::Named( "IM" ) = meanBlocks );

}

Rcpp::List meanByBlocks( const Array & M, const IVector & clu,  const IVector & nClu, const IVector & n, const std::string diagonal,
                         const std::string & sBorders, const Rcpp::Nullable<Array> & bordersMatLower, const Rcpp::Nullable<Array> & bordersMatUpper,
                         const Rcpp::Nullable<DMatrix> & bordersSeperateLower, const Rcpp::Nullable<DMatrix> & bordersSeperateUpper, const bool & addOne, const double & eps)
{
	const int dimensions = Rcpp::sum( nClu );
		
    Diagonale dDiag = getDiagonale( diagonal );
    const BorderType eBorders = getBorderType( sBorders );
    Borders<Array> bordersMeanstMat;
    Borders<DMatrix> bordersSeperate;
	
	IVector unitMode(clu.size());
	int ind=0;
	for(int i = 0; i<clu.size();++i){
		if(i==n.at(ind)){
			++ind;
		}
		unitMode.at(i)=ind;
	}
	
//    if( useBorders ) {
//        checkInputBorders( dDiag, bordersMatLower, bordersMatUpper, bordersSeperateLower, bordersSeperateUpper );
//        if( dDiag == Diagonale::Seperate ) {
//            bordersSeperate = Borders<DMatrix>( Rcpp::as<DMatrix>( bordersSeperateLower ), Rcpp::as<DMatrix>( bordersSeperateUpper ) );
//        }
//        bordersMeanstMat = Borders<Array>( Rcpp::as<Array>( bordersMatLower ), Rcpp::as<Array>( bordersMatUpper ) );
//    }
//    else {
//        if( dDiag == Diagonale::Seperate ) {
//            bordersSeperate = Borders<DMatrix>( dimensions, M.n_slices );
//        }
//        bordersMeanstMat = Borders<Array>( dimensions, dimensions, M.n_slices );
//    }

    if( eBorders != BorderType::None ) {
        checkInputBorders( dDiag, bordersMatLower, bordersMatUpper, bordersSeperateLower, bordersSeperateUpper );
        if( dDiag == Diagonale::Seperate ) {
            bordersSeperate = Borders<DMatrix>( Rcpp::as<DMatrix>( bordersSeperateLower ), Rcpp::as<DMatrix>( bordersSeperateUpper ) );
        }
        bordersMeanstMat = Borders<Array>( Rcpp::as<Array>( bordersMatLower ), Rcpp::as<Array>( bordersMatUpper ) );
    }
    else {
        if( dDiag == Diagonale::Seperate ) {
            bordersSeperate = Borders<DMatrix>( dimensions, M.n_slices );
        }
        bordersMeanstMat = Borders<Array>( dimensions, dimensions, M.n_slices );
    }

    DMatrix pSeparate;
    Array aRes;
	const size_t S = static_cast<size_t>( n.size() );
	Array superBlockMeans(S, S, M.n_slices, arma::fill::zeros );
	DMatrix superBlockDiagMeans( S, M.n_slices, arma::fill::zeros );

	superblockMeansFun(M, n, dDiag, unitMode, superBlockMeans, superBlockDiagMeans);
//    const DMatrix MEANS = relationsMeans( M, n );
    meansByBlocks( M, aRes, clu, nClu, pSeparate, superBlockMeans, superBlockDiagMeans, n, eBorders, bordersMeanstMat, bordersSeperate, dDiag, addOne, eps);

    if( !pSeparate.is_empty() ) {
        return Rcpp::List::create( Rcpp::Named( "meansByBlocs" ) = aRes, Rcpp::Named( "meansByCluDiag" ) = pSeparate );
    }
    return Rcpp::List::create( Rcpp::Named( "meansByBlocs" ) = aRes );

}

double critFunction( const Array & M, const IVector & clu, const Array & weights, const DVector & uWeights, const int dimensions, const IVector & n, const double weightClusterSize, const std::string & diagonal,
                     const std::string & sBorders, const Rcpp::Nullable<Array> & bordersMatLower, const Rcpp::Nullable<Array> & bordersMatUpper,
                     const Rcpp::Nullable<DMatrix> & bordersSeperateLower, const Rcpp::Nullable<DMatrix> & bordersSeperateUpper, const bool & addOne, const double & eps)
{
    Diagonale dDiagonale = getDiagonale( diagonal );
    const BorderType eBorders = getBorderType( sBorders );
    Borders<Array> bordersMeanstMat;
    Borders<DMatrix> bordersSeperate;
	
	const int K = dimensions;
	IVector cumN = Rcpp::cumsum( n );
	IVector countGroups( K ); 
	IVector cluMode( K );
	IVector nClu(n.size(),0);
	int ind = 0;
	for( int i = 0; i < static_cast<int>( clu.size() ); ++i ) {
		++countGroups.at(clu.at(i));
		if(i==cumN.at(ind)){
			++ind;
		}
		cluMode.at(clu.at(i))=ind;
		//Rcpp::Rcout << "i = " << i << ", cluMode(" << clu.at(i) << ") = "<< ind << "\n";
	}
	
	for(int i = 0; i < static_cast<int>( K ); ++i ) {
		++nClu.at(cluMode.at(i));
	}
	
	// Rcpp::Rcout << "nClu:\n" << nClu << "\n";
	// Rcpp::Rcout << "cluMode:\n" << cluMode << "\n";
	
	// for(int i = 0;i<n.size();++i){
		// for(int j = 0; j<nClu.at(i);++j){
			// cluMode.at(ind)=i;
			// ++ind;
			// Rcpp::Rcout << "cluMode(" << ind <<") = "<< i << "\n";
		// }
	// }
	DVector logProbGroups( K );
	for(int i = 0; i< K;++i){
		logProbGroups.at(i)=log(static_cast<double>(countGroups.at(i))/n.at(cluMode.at(i)));
	}
	
	// Rcpp::Rcout << "logProbGroups: " << logProbGroups <<"\n";
	// Rcpp::Rcout << "clu.size(): " << clu.size() <<"\n";
	IVector unitMode(clu.size());
	// Rcpp::Rcout << "unitMode: " << unitMode <<"\n";
	// Rcpp::Rcout << "n: " << n <<"\n";
	ind=0;
	int nNext = n.at(ind);
	for(int i = 0; i<clu.size();++i){
		// Rcpp::Rcout << "i: " << i <<"\n";
		if(i==nNext){
			++ind;
			nNext = nNext + n.at(ind);
		}
		unitMode.at(i)=ind;
	}	
	// Rcpp::Rcout << "unitMode: " << unitMode <<"\n";
	// Rcpp::Rcout << "OK1\n";
//    if( useBorders ) {
//        checkInputBorders( dDiagonale, bordersMatLower, bordersMatUpper, bordersSeperateLower, bordersSeperateUpper );
//        if( dDiagonale == Diagonale::Seperate ) {
//            bordersSeperate = Borders<DMatrix>( Rcpp::as<DMatrix>( bordersSeperateLower ), Rcpp::as<DMatrix>( bordersSeperateUpper ) );
//        }
//        bordersMeanstMat = Borders<Array>( Rcpp::as<Array>( bordersMatLower ), Rcpp::as<Array>( bordersMatUpper ) );
//    }
//    else {
//        if( dDiagonale == Diagonale::Seperate ) {
//            bordersSeperate = Borders<DMatrix>( dimensions, M.n_slices );
//        }
//        bordersMeanstMat = Borders<Array>( dimensions, dimensions, M.n_slices );
//    }

    if( eBorders != BorderType::None ) {
        checkInputBorders( dDiagonale, bordersMatLower, bordersMatUpper, bordersSeperateLower, bordersSeperateUpper );
        if( dDiagonale == Diagonale::Seperate ) {
            bordersSeperate = Borders<DMatrix>( Rcpp::as<DMatrix>( bordersSeperateLower ), Rcpp::as<DMatrix>( bordersSeperateUpper ) );
        }
        bordersMeanstMat = Borders<Array>( Rcpp::as<Array>( bordersMatLower ), Rcpp::as<Array>( bordersMatUpper ) );
    }
    else {
        if( dDiagonale == Diagonale::Seperate ) {
            bordersSeperate = Borders<DMatrix>( dimensions, M.n_slices );
        }
        bordersMeanstMat = Borders<Array>( dimensions, dimensions, M.n_slices );
    }

	// Rcpp::Rcout << "OK2\n";
	
    DMatrix pSeparate;
    Array aRes;
//  const DMatrix MEANS = relationsMeans( M, n );
	const size_t S = static_cast<size_t>( n.size() );
	Array superBlockMeans(S, S, M.n_slices, arma::fill::zeros );
	DMatrix superBlockDiagMeans( S, M.n_slices, arma::fill::zeros );
	superblockMeansFun(M, n, dDiagonale, unitMode, superBlockMeans, superBlockDiagMeans);
	// Rcpp::Rcout << "OK3\n";
	
    meansByBlocks( M, aRes, clu, nClu, pSeparate, superBlockMeans, superBlockDiagMeans, n, eBorders, bordersMeanstMat, bordersSeperate, dDiagonale, addOne, eps);
	// Rcpp::Rcout << "OK4\n";

	// Rcpp::Rcout << "clu\n" << clu <<"\n";		
	// Rcpp::Rcout << "M\n" << M <<"\n";		
	// Rcpp::Rcout << "weights\n" << weights <<"\n";	
	// Rcpp::Rcout << "meanBlocks\n" << aRes <<"\n";		
	// Rcpp::Rcout << "pSeparate\n" << pSeparate <<"\n";		
	// Rcpp::Rcout << "dDiag\n" << dDiagonale <<"\n";
	// Rcpp::Rcout << "logProbGroups\n" << logProbGroups <<"\n";
	// Rcpp::Rcout << "weightClusterSize\n" << weightClusterSize <<"\n";
	
    return criterialFunction( M, clu, weights, uWeights, aRes, pSeparate, dDiagonale, logProbGroups, weightClusterSize);
}

void meansByBlocks( const Array & M, Array & res, const IVector & clu, const IVector & nClu, DMatrix & p_pSepare, const Array & superBlockMeans, const DMatrix & superBlockDiagMeans, const IVector & n, BorderType p_type,
                    const Borders<Array> & p_btBorders, const Borders<DMatrix> & p_btBordersSeperate, const Diagonale sDiagonal, const bool & addOne, const double & eps)
{
	// Rcpp::Rcout << "nClu = " << nClu << "\n";
	const int dimensions = Rcpp::sum(nClu);
	double dVal;
	// Rcpp::Rcout << "dimensions = " << dimensions << "\n";

    if( res.is_empty() ) {
        res = Array( dimensions, dimensions, M.n_slices, arma::fill::zeros );
    }
    else {
        res.fill( 0 );
    }
    Array S = res;
    Array N = res;
    DMatrix mDiagonalRes;
    DMatrix mSseprateDiagonal;
    DMatrix mNseprateDiagonal;

    if( sDiagonal == Diagonale::Seperate ) {
        mDiagonalRes = DMatrix( dimensions, M.n_slices, arma::fill::zeros );
        mSseprateDiagonal = mDiagonalRes;
        mNseprateDiagonal = mDiagonalRes;
    }
	// Rcpp::Rcout << "MB1\n";
	
    for( size_t r = 0; r < M.n_slices; ++r ) {
        for( size_t i = 0; i < M.n_rows; ++i ) {
            for( size_t j = 0; j < M.n_cols; ++j ) {
                if( sDiagonal == Diagonale::Ignore && i == j ) { // Ignore diagonal
                    continue;
                }
                else if( sDiagonal == Diagonale::Seperate && i == j ) { // Calculate diagonal seperately and ignore it in return Array
                    mSseprateDiagonal( clu.at( i ), r ) += M( i, j, r );
                    mNseprateDiagonal( clu.at( i ), r ) += 1;
                    continue;
                }
                S( clu.at( i ), clu.at( j ), r ) += M( i, j, r );
                N( clu.at( i ), clu.at( j ), r ) += 1;
            }
        }
    }

	// Rcpp::Rcout << "Within function meansByBlocks\n";
	// Rcpp::Rcout << "sums (S) \n";
	// Rcpp::Rcout << S <<"\n";
	// Rcpp::Rcout << "\n";
	
	// Rcpp::Rcout << "number of cells (N) \n";
	// Rcpp::Rcout << N <<"\n";
	// Rcpp::Rcout << "\n";
	
//    for( size_t r = 0; r < res.n_slices; ++r ) {
//        int s = 0;
//        int sCount = 0;
//        for( size_t i = 0; i < res.n_rows; ++i ) {
//            if( sDiagonal == Diagonale::Seperate ) { // save diagonal values into Matrix[ dimensions, r ]
//                double dVal = double( mSseprateDiagonal( i, r ) ) / mNseprateDiagonal( i, r );
//                if( dVal < p_btBordersSeperate.getLowerAt( i, r ) ) mDiagonalRes.at( i, r ) = p_btBordersSeperate.getLowerAt( i, r );
//                if( dVal > p_btBordersSeperate.getUpperAt( i, r ) ) mDiagonalRes.at( i, r ) = p_btBordersSeperate.getUpperAt( i, r );
//                else mDiagonalRes( i, r ) = dVal;
//            }
//            if( sCount >= n.at( s ) ) {
//                ++s;
//                sCount = 0;
//            }
//            ++sCount;
//            for( size_t j = 0; j < res.n_cols; ++j ) {
//                double dVal( S( i, j, r ) );
//                if( i == j && N( i, j, r ) == 0 && ( ( sDiagonal == Diagonale::Ignore ) || ( sDiagonal == Diagonale::Seperate ) ) ) { // If value of the block is 0 and we ignored diagonal values, set value of the block to mean (M[ , , r ] )
//                    res(i, j, r ) = p_mMeans.at( s, r );
//                }
//                else {
//                    res( i, j, r ) = dVal / N( i, j, r );
//                }

//                if( res.at( i, j, r ) < p_btBorders.getLowerAt( i, j, r ) ) res.at( i, j, r ) = p_btBorders.getLowerAt( i, j, r );
//                else if( res.at( i, j, r ) > p_btBorders.getUpperAt( i, j, r ) ) res.at( i, j, r ) = p_btBorders.getUpperAt( i, j, r );
//            }
//        }
//    }
// superBlockMeans
//superBlockDiagMeans
	// Rcpp::Rcout << "Within function meansByBlocks - test determining of sets\n";

    for( size_t r = 0; r < res.n_slices; ++r ) {
        int s1 = 0;
        int sCount1 = 0;
        for( size_t i = 0; i < res.n_rows; ++i ) {
            if( sCount1 >= nClu.at( s1 ) ) {
                ++s1;
                sCount1 = 0;
            }
            ++sCount1;
            if( sDiagonal == Diagonale::Seperate ) { // save diagonal values into Matrix[ dimensions, r ]
				if(addOne){
					dVal = double( mSseprateDiagonal( i, r ) +  superBlockDiagMeans.at(s1,r)) / (mNseprateDiagonal( i, r ) + 1);
				} else {
					dVal = double( mSseprateDiagonal( i, r )) / mNseprateDiagonal( i, r );
					if (dVal < eps){
						dVal = eps;
					} else if(dVal > (1-eps)){
						dVal = 1-eps;
					}
				}

                mDiagonalRes( i, r ) = dVal;
            }
			
			int s2 = 0;
			int sCount2 = 0;	
            for( size_t j = 0; j < res.n_cols; ++j ) {
				if( sCount2 >= nClu.at( s2 ) ) {
					++s2;
					sCount2 = 0;
				}
				++sCount2;

                double dVal( S( i, j, r ) );
                if( i == j && N( i, j, r ) == 0 && ( ( sDiagonal == Diagonale::Ignore ) || ( sDiagonal == Diagonale::Seperate ) ) ) { // If value of the block is 0 and we ignored diagonal values, set value of the block to mean (M[ , , r ] )
                    res(i, j, r ) = superBlockMeans.at( s1, s1, r );
                }
                else {
					if(addOne){
						res( i, j, r ) = (dVal + superBlockMeans(s1,s2,r))/ (N( i, j, r )+1);
					} else {
						res( i, j, r ) = dVal/N( i, j, r );
						if (res( i, j, r ) < eps){
							res( i, j, r ) = eps;
						} else if(res( i, j, r ) > (1-eps)){
							res( i, j, r ) = 1-eps;
						}
					}
                }
				// Rcpp::Rcout << "clu1 = " << i <<", clu2 = " << j << ", set1 = " << s1 << ", set2 = " << s2<< "\n";
            }
        }
    }


	// Rcpp::Rcout << "Within function meansByBlocks\n";
	// Rcpp::Rcout << "meanBlocks - before applying borders \n";
	// Rcpp::Rcout << res <<"\n";
	// Rcpp::Rcout << "\n";	

    if( p_type != BorderType::None ) {
        for( size_t r = 0; r < res.n_slices; ++r ) {
            for( size_t i = 0; i < res.n_rows; ++i ) {
                if( sDiagonal == Diagonale::Seperate ) { // save diagonal values into Matrix[ dimensions, r ]
                    double dVal = double( mSseprateDiagonal( i, r ) ) / mNseprateDiagonal( i, r );
                    if( p_type == BorderType::Inside ) {
                        if( dVal < p_btBordersSeperate.getLowerAt( i, r ) ) mDiagonalRes.at( i, r ) = p_btBordersSeperate.getLowerAt( i, r );
                        if( dVal > p_btBordersSeperate.getUpperAt( i, r ) ) mDiagonalRes.at( i, r ) = p_btBordersSeperate.getUpperAt( i, r );
                        else mDiagonalRes( i, r ) = dVal;
                    }
                    else if( p_type == BorderType::Outside ) {
                        if( dVal > p_btBordersSeperate.getLowerAt( i, r ) && dVal < p_btBordersSeperate.getUpperAt( i, r ) ) {
                            double dDiffLower = std::abs( p_btBordersSeperate.getLowerAt( i, r ) - dVal );
                            double dDiffUpper = std::abs( p_btBordersSeperate.getUpperAt( i, r ) - dVal );
                            if( dDiffLower - dDiffUpper < 0 ) mDiagonalRes.at( i, r ) = p_btBordersSeperate.getLowerAt( i, r );
                            else mDiagonalRes.at( i, r ) = p_btBordersSeperate.getUpperAt( i, r );
                        }
                    }
                }
                for( size_t j = 0; j < res.n_cols; ++j ) {
                    double dVal = res.at( i, j, r );
                    if( p_type == BorderType::Inside ) {
                        if( dVal < p_btBorders.getLowerAt( i, j, r ) ) res.at( i, j, r ) = p_btBorders.getLowerAt( i, j, r );
                        else if( dVal > p_btBorders.getUpperAt( i, j, r ) ) res.at( i, j, r ) = p_btBorders.getUpperAt( i, j, r );
                    }
                    else if( p_type == BorderType::Outside ) {
                        if( dVal > p_btBorders.getLowerAt( i, j, r ) && dVal < p_btBorders.getUpperAt( i, j, r ) ) {
                            double dDiffLower = std::abs( p_btBorders.getLowerAt( i, j, r ) - dVal );
                            double dDiffUpper = std::abs( p_btBorders.getUpperAt( i, j, r ) - dVal );
                            if( dDiffLower - dDiffUpper < 0 ) res.at( i, j, r ) = p_btBorders.getLowerAt( i, j, r );
                            else res.at( i, j, r ) = p_btBorders.getUpperAt( i, j, r );
                        }
                    }
                }
            }
        }
    }
	
	// Rcpp::Rcout << "Within function meansByBlocks\n";
	// Rcpp::Rcout << "meanBlocks - after applying borders \n";
	// Rcpp::Rcout << res <<"\n";
	// Rcpp::Rcout << "\n";	

    if( sDiagonal == Diagonale::Seperate ) { // Save seperate digaonal values to input parameter
        p_pSepare = std::move( mDiagonalRes );
    }
}

double criterialFunction( const Array & M, const IVector & clu, const Array & weights, const DVector & uWeights, const Array & meansMat, const DMatrix & p_mSeparate, const Diagonale p_diagonale, const DVector & logProbGroups, const double weightClusterSize)
{
	// Rcpp::Rcout << "CF0\n";	
    double dRet = 0;
	// double dTmp = 0;
    for( size_t i = 0; i < M.n_rows; ++i ) {
		// Rcpp::Rcout << "i = " << i <<"\n";	
		dRet += -uWeights.at(i)*logProbGroups.at(clu.at( i ))*weightClusterSize;
        for( size_t j = 0; j < M.n_cols; ++j ) {
			// if(i==0) Rcpp::Rcout << "j = " << j <<"\n";	
            for( size_t r = 0; r < M.n_slices; ++r ) {
				// if(i==0 && j==0) Rcpp::Rcout << "r = " << r <<"\n";	
				if(weights( i, j, r )<=0) continue;
                if( p_diagonale == Diagonale::Ignore && i == j ) { // ignore diagonal
                    continue;
                }
                else if( p_diagonale == Diagonale::Seperate && i == j ) {
                    double dAvg = p_mSeparate( clu.at( i ), r );                   
					// dTmp = -weights( i, j, r ) * (M( i, j, r )* log(dAvg) +  (1 - M( i, j, r ))*log(1 - dAvg));
					dRet += -weights( i, j, r ) * (M( i, j, r )* log(dAvg) +  (1.0 - M( i, j, r ))*log(1.0 - dAvg));
                }
                else {
                    double dAvg = meansMat( clu.at( i ), clu.at( j ), r );
                    // dTmp = -weights( i, j, r ) * (M( i, j, r )* log(dAvg) +  (1 - M( i, j, r ))*log(1 - dAvg));
					dRet += -weights( i, j, r ) * (M( i, j, r )* log(dAvg) +  (1.0 - M( i, j, r ))*log(1.0 - dAvg));
                }
				// Rcpp::Rcout << "LL(" << i << ", " << j << ", " << r << ") = " << dTmp <<"\n";	
            }
        }
    }

    return  dRet;
}

void setGroups( const Array & M, IVector & clu, const Array & weights, const DVector & uWeights, const Array & meansMat, const IVector & nClu, const IVector & n, const DMatrix & p_mSeparate, const Diagonale p_diagonale, const double weightClusterSize, const IVector & cluMode,IVector & countGroups, DVector & logProbGroups)
{
    IVector borders = Rcpp::cumsum( nClu );
    double eMin;
    double eTmp;
	int nkMin = 1;
    int K = Rcpp::sum( nClu );
    size_t kMin( 0 );
    DVector eVec( clu.size() );
    IVector countGroupsNew( K );
    DVector e( K );
    for( unsigned int i = 0; i < static_cast<unsigned int>( clu.size() ); ++i ) {
        eMin = DBL_MAX;
		nkMin=0;
        int group( clu.at( i ) );
        size_t iBelongsTo( belongsTo( group, borders ) );
        size_t k = 0;
        if( iBelongsTo ) {
            k = borders.at( iBelongsTo - 1 );
        }
        for( ; k < static_cast<unsigned int>( borders.at( iBelongsTo ) ); ++k ) {
            eTmp = 0;
            for( unsigned int j = 0; j < static_cast<unsigned int>( clu.size() ); ++j ) {
                size_t cluJ = static_cast<size_t>( clu.at( j ) );
                if( i == j ){
					if( p_diagonale == Diagonale::Ignore) { 
						continue; 
					} else if( p_diagonale == Diagonale::Seperate ) {
						for( unsigned int r = 0; r < M.n_slices; ++r ) {
							if(weights( i, j, r )<=0) continue;
                            double dAvg = p_mSeparate.at( k, r );
                            eTmp += - weights( i, j, r ) * (M( i, j, r )* log(dAvg) +  (1 - M( i, j, r ))*log(1 - dAvg));
						}
					} else if( p_diagonale == Diagonale::Same ) {
						for( unsigned int r = 0; r < M.n_slices; ++r ) {
							if(weights( i, j, r )<=0) continue;
                            double dAvg = meansMat( k, k, r );
                            eTmp += - weights( i, j, r ) * (M( i, j, r )* log(dAvg) +  (1 - M( i, j, r ))*log(1 - dAvg));
						}
					}
				} else {
                    for( unsigned int r = 0; r < M.n_slices; ++r ) {
						if(weights( i, j, r )>0){
							double dAvg = meansMat( k, cluJ, r );
							eTmp += - weights( i, j, r ) *(M( i, j, r )* log(dAvg) +  (1 - M( i, j, r ))*log(1 - dAvg));	
						}
						if(weights( j, i, r )>0){
							double dAvg = meansMat( cluJ, k, r );
							eTmp += -weights( j, i, r ) * (M( j, i, r )* log(dAvg) +  (1 - M( j, i, r ))*log(1 - dAvg));
						}
						
					
                    }
                }

            }
			eTmp += -uWeights.at(i)*logProbGroups.at(k)*weightClusterSize;
			// Rcpp::Rcout << "i = " << i <<", k = " << k << ", eTmp= " << eTmp << "\n";
            if ( eTmp < eMin ){
                kMin = k;
                eMin = eTmp;
				nkMin = 1;
            } else if (eTmp == eMin){
				++nkMin;
				if(R::runif(0,1)< 1.0/nkMin){
					kMin = k;
				}
			}

        }
        clu.at( i ) = static_cast<int>( kMin );
        eVec.at( i ) = eMin;
        ++countGroupsNew.at( kMin );
    }

    IVector nBorders( Rcpp::cumsum( n ) );


    int iBegin, iEnd, k;
    for( size_t i = 0; i < static_cast<size_t>( nBorders.size() ); ++i ){
        if( !i ) {
            iBegin = 0;
            k = 0;
        }
        else {
            iBegin = nBorders.at( i - 1 );
            k = borders.at( i - 1 );
        }
        iEnd = nBorders.at( i );
        K = borders.at( i );
        for( ; k < K; ++k ) {
            if( countGroupsNew.at(k)==0) {
                size_t g = std::distance( eVec.begin(), std::max_element( eVec.begin() + iBegin, eVec.begin() + iEnd ) );
//                ce ma countgroups[i] samo 1 skupino, zberem naslednji max iz drugih skupin - torej ce je k 1 - 3 in ima skupina 1 samo 1 skupino izberem max med skupinama 2 in 3
                if( countGroupsNew.at( clu.at( g ) ) < 2 ) {
                    k--;
                    eVec.at( g ) = -1;
                    continue;
                }
				--countGroupsNew.at(clu.at( g ));
				++countGroupsNew.at(k);
                clu.at( g ) = k;
                eVec.at( g ) = 0;
                k = i == 0 ? 0 : borders.size() > 1 ? borders.at( i - 1 ) : 0;
                --k;
            }
        }
    }
	
	for(int i = 0; i< K;++i){
		countGroups.at(i)=countGroupsNew.at(i);
	}
	
	for(int i = 0; i< K;++i){
		logProbGroups.at(i)=log(static_cast<double>(countGroups.at(i))/n.at(cluMode.at(i)));
	}
	
	

}

unsigned int belongsTo( const int & group, const IVector & borders )
{
    for( size_t i = 0; i < static_cast<size_t>( borders.size() ); ++i ) {
        if( group < borders.at( i ) ) {
            return i;
        }
    }
    return static_cast<size_t>( Rcpp::sum( borders ) );
}


DMatrix relationsMeans( const Array & M, const IVector & n )
{
    const size_t S = static_cast<size_t>( n.size() );
    DMatrix mMeans( S, M.n_slices, arma::fill::zeros );
    IVector cumN = Rcpp::cumsum( n );
    cumN.push_front( 0 );
    for( size_t s = 0; s < S; ++s ) {
        for( size_t r = 0; r < M.n_slices; ++r ) {
            Array subArray = M.subcube( cumN.at( s ), cumN.at( s ), r, cumN.at( s + 1 ) - 1, cumN.at( s + 1 ) - 1, r );
            mMeans.at( s, r ) = meanMatrix( subArray.slice( 0 ) );
        }
    }

    return mMeans;
}


double meanMatrix( const DMatrix & p_matrix )
{
    size_t sElementsN = 0;
    double sum = 0;
    for( size_t i = 0; i < static_cast<size_t>( p_matrix.n_rows ); ++i ) {
        for( size_t j = 0; j < static_cast<size_t>( p_matrix.n_cols ); ++j ) {
            if( i == j ) continue;
            sum += p_matrix.at( i, j );
            ++sElementsN;
        }
    }
    return sum/ sElementsN;
}


void superblockMeansFun( const Array & M, const IVector & n, const Diagonale p_diagonale, const IVector & unitMode, Array & superBlockMeans, DMatrix & superBlockDiagMeans)
{
    const size_t S = static_cast<size_t>( n.size() );
    IVector cumN = Rcpp::cumsum( n );
    cumN.push_front( 0 );
	size_t sElementsN;
	double sum;
	size_t sDiagElementsN;
	double diagSum;	
    for( size_t s1 = 0; s1 < S; ++s1 ) {
		for( size_t s2 = 0; s2 < S; ++s2 ) {
			for( size_t r = 0; r < M.n_slices; ++r ) {
				Array subArray = M.subcube( cumN.at( s1 ), cumN.at( s2 ), r, cumN.at( s1 + 1 ) - 1, cumN.at( s2 + 1 ) - 1, r );
				DMatrix p_matrix = subArray.slice(0);
				sElementsN = 0;
				sum = 0;
				sDiagElementsN = 0;
				diagSum = 0;

				for( size_t i = 0; i < static_cast<size_t>( p_matrix.n_rows ); ++i ) {
					for( size_t j = 0; j < static_cast<size_t>( p_matrix.n_cols ); ++j ) {
						if( i == j && s1==s2 ){
							if( p_diagonale == Diagonale::Ignore) { 
							continue; 
							}
							else if( p_diagonale == Diagonale::Seperate ) {
								diagSum += p_matrix.at( i, j );
								++sDiagElementsN;					
							}
							else if( p_diagonale == Diagonale::Same ) {
								sum += p_matrix.at( i, j );
								++sElementsN;
							}
						} else{
							sum += p_matrix.at( i, j );
							++sElementsN;
						}
					}
				}
				superBlockMeans(s1,s2,r)= sum/ sElementsN;

				if(s1==s2 ){
					superBlockDiagMeans.at(s1,r)=diagSum/sDiagElementsN;
				}
					
			}
		}
    }
}



Diagonale getDiagonale( const std::string & p_sDiagonal )
{
    Diagonale dRet = Diagonale::Ignore;
    std::string sTmp = p_sDiagonal;
    std::transform( sTmp.begin(), sTmp.end(), sTmp.begin(), ::tolower );
    if( sTmp == "ignore" );
    else if( sTmp == "same" ) {
        dRet = Diagonale::Same;
    }
    else if( sTmp == "seperate" ) {
        dRet = Diagonale::Seperate;
    }
    else Rcpp::stop( "Unknow diagonal parameter\nOptions are: [ same, ignore, seperate ]\n" );

    return dRet;
}

BorderType getBorderType( const std::string & p_sBorder )
{
    BorderType eRet = BorderType::None;

    std::string sTmp = p_sBorder;
    std::transform( sTmp.begin(), sTmp.end(), sTmp.begin(), tolower );
    if( sTmp == "none" );
    else if( sTmp == "inside" ) {
        eRet = BorderType::Inside;
    }
    else if( sTmp == "outside" ) {
        eRet = BorderType::Outside;
    }
    else {
        Rcpp::stop( "Unknown border type\nOptions are: [ none, inside, seperate ]\n" );
    }

    return eRet;
}

void checkInputBorders( const Diagonale & p_diagonale, const Rcpp::Nullable<Array> & bordersMatLower, const Rcpp::Nullable<Array> & bordersMatUpper,
                    const Rcpp::Nullable<DMatrix> & bordersSeperateLower, const Rcpp::Nullable<DMatrix> & bordersSeperateUpper )
{
    if( bordersMatLower.isNull() ) Rcpp::stop( "Invalid argument: bordersMatLower is null" );
    if( bordersMatUpper.isNull() ) Rcpp::stop( "Invalid argument: bordersMatUpper is null" );
    if( p_diagonale == Diagonale::Seperate ) {
        if( bordersSeperateLower.isNull() ) Rcpp::stop( "Invalid argument: bordersSeperateLower is null" );
        if( bordersSeperateUpper.isNull() ) Rcpp::stop( "Invalid argument: bordersSeperateUpper is null" );
    }
}

std::ostream & operator << ( std::ostream & p_stream, const Diagonale p_diag )
{
    switch ( p_diag ) {
        case Diagonale::Ignore:
            p_stream << "Ignore";
        break;
        case Diagonale::Same:
            p_stream << "Same";
        break;
        case Diagonale::Seperate:
            p_stream << "Seperate";
        break;
    }

    return p_stream;
}
