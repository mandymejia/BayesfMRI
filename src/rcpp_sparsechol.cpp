using Eigen::Map;
using Eigen::SparseMatrix;
using Eigen::LLT;
const SparseMatrix<double> SS(as<SparseMatrix<double> >(Sigma));

typedef Eigen::SimplicialLLT<SparseMatrix<double> > SpChol;

const SpChol Ch(SS);
return wrap(Ch.matrixL());
