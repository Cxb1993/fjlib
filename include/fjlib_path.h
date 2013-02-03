#ifndef TFJLibDirectoryH
#define TFJLibDirectoryH

#include "fjlib_string.h"

#ifdef __GNUC__
#include <unistd.h>
#include <sys/stat.h>
#else
#include <io.h>
#include <direct.h>
#endif

namespace fjlib {

#define MAX_STR_LEN 256

// http://flinflon.brandonu.ca/Dueck/SystemsProgramming/MultiShell/Directory.htm
class TFJPath {
private:
	void		clear()
	{
		memset(_drive,0,sizeof(_drive));
		memset(_dir,0,sizeof(_dir));
		memset(_fname,0,sizeof(_fname));
		memset(_ext,0,sizeof(_ext));
	}
protected:
    char		*_drive,*_dir,*_fname,*_ext;
public:
	TFJPath() 
	{
		_drive=new char[MAX_STR_LEN];
		_dir=new char[MAX_STR_LEN];
		_fname=new char[MAX_STR_LEN];
		_ext=new char[MAX_STR_LEN];
	}
	~TFJPath()
	{
		delete[] _ext;
		delete[] _fname;
		delete[] _dir;
		delete[] _drive;
	}
	#ifndef __GNUC__	// doesn't work on linux
	bool		splitpath(const str_t& path)
	{
		if (path=="") return false;
		clear();
		_splitpath(path.c_str(),_drive,_dir,_fname,_ext);
		return true;
	}
	#endif
	const str_t	drive()	const	{ return str_t(_drive); }
	const str_t	dir() const		{ return str_t(_dir); }
	const str_t	filename() const	{ return str_t(_fname); }
	const str_t	ext() const		{ return str_t(_ext); }

	static bool exists(const str_t& path)
	{
		#ifdef __GNUC__ 
			return (::access(path.c_str(),0)!=-1);
		#else
			return (_access(path.c_str(),0)!=-1);
		#endif
	}
	static bool chdir(const str_t& dir)
	{
		#ifdef __GNUC__
			return (::chdir(dir.c_str())==0);
		#else
			return (_chdir(dir.c_str())==0);
		#endif
	}
	static bool mkdir(const str_t& dir)
	{
		#ifdef __GNUC__
			return (::mkdir(dir.c_str(),0755)==0);
		#else
			return (_mkdir(dir.c_str())==0);
		#endif
	}
	static bool rmdir(const str_t& dir)
	{
		#ifdef __GNUC__
			return ::rmdir(dir.c_str());
		#else
			return  _rmdir(dir.c_str());
		#endif
	}
	static
	const str_t	pwd()
	{
		char dirname[MAX_STR_LEN];
		#ifdef __GNUC__
			char* rc=::getcwd(dirname, sizeof(dirname));
		#else
			char* rc=_getcwd(dirname, sizeof(dirname));
		#endif
		if (rc==NULL) return str_t("");
		else return str_t(rc);
	}
};


}

#endif
