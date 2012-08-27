#ifndef __param_tags__
#define __param_tags__
/*!
  \brief Defines standard string tags for various cosmological parameters
*/
#include<string>
namespace param_tags {

  //paramcodes must be the unambiguous parameter name, where
  // the tags below are used for human readable output and so
  // don't have to be unique
  //The order should be that we output prob surfaces in
  enum paramcodes { omegam=0, omegade=1, w0=2, wa=3, alpha=100, beta=101, 
		    scriptm=102, unknown = 10000 }; //!< Codes for various parameters
  enum paramtype { unset, nuisance, cosmological }; //!< What type of parameter

  typedef std::pair< paramcodes, paramtype > paramspec;

  const std::string omtag = "\\Omega_{m}"; //!< \f$\Omega_m\f$
  const std::string oltag = "\\Omega_{\\Lambda}"; //!< \f$\Omega_{\Lambda}\f$
  const std::string odetag = "\\Omega_{DE}"; //!< \f$\Omega_{DE}\f$
  const std::string wtag = "w"; //!< \f$w\f$
  const std::string w0tag = "w_0"; //!< \f$w_0\f$
  const std::string watag = "w_a"; //!< \f$w_a\f$
  const std::string alphatag = "\\alpha"; //!< \f$\alpha\f$
  const std::string betatag = "\\beta"; //!< \f$\beta\f$
  const std::string scriptmtag = "\\mathcal{M}"; //!< \f${\mathcal M}\f$
}
#endif
