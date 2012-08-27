#ifndef __cosfitterexcept__
#define __cosfitterexcept__

#include <string>
#include <iostream>

/*! 
  \brief Exception class for simple_cosfitter 
*/

class CosFitterExcept {

  void init(const std::string& errclass,const std::string& errmethod,
	    const std::string& errstr,int err);  //!< Internal initialization
  bool classset; //!< Is errclass set
  bool methodset; //!< Is errmethod set
  bool strset; //!< Is errstr set
  bool errset; //!< Is errnum set

 public:
  std::string errclass;          //!< Class throwing the exception
  std::string errmethod;         //!< Method throwing the exception
  std::string errstr;            //!< Error string (user consumption)
  int errnum;               //!< Error code (class or method specific)

  // Constructors

  CosFitterExcept(); //!< Basic constructor
  explicit CosFitterExcept(const std::string errstr); //!< Just with errstring
  CosFitterExcept(const std::string errstr,int err); //!< Errstring and number
  CosFitterExcept(const std::string errclass,const std::string errmethod,
		  const std::string errstr); //!< Class, method, error string
  CosFitterExcept(const std::string errclass,const std::string errmethod,
		  const std::string errstr,int err); //!< Class, method, error string, and number

  friend std::ostream& operator<<(std::ostream& os,const CosFitterExcept& ex);//!< Output operator for CosFitterExcept

  
};

#endif
