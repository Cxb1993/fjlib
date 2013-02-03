#ifndef FJLIB_LOG
#define FJLIB_LOG

#include "fjlib_string.h"
#include <fstream>

// 12.30	created

namespace fjlib {

/*!
//	TFJLog is used for logging message to the memory
//	and the message can be saved to a file whenever needed.
//	Usage:
//	TFJLog log;
//	log << "first line\n";
//	log.SaveLog();
//	
!*/
class TFJLog {
private:
	str_t		_filename;
	std::ostream*
				_os;
	bool		_save_on_destroy;
protected:
	std::ostringstream
				_cache;
public:
	TFJLog(): _filename(""), _os(NULL), _save_on_destroy(true) {}
	TFJLog(const str_t& filename, std::ostream* os=NULL): _filename(filename), _os(os), _save_on_destroy(true)  {}
	~TFJLog() { if (_save_on_destroy) SaveLog(false); }
	void		SetFilename(const str_t& fn) { _filename=fn; }
	void		UseExtraStream(std::ostream* os) { _os=os; }
	void		NoExtraStream() { _os=NULL; }
	void		SetSaveOnDestroy(bool save) { _save_on_destroy=save; }
	// add log
/*
	void		AddLog(const str_t& log_string)
	{
		_cache << log_string;
		if (_os!=NULL) (*_os) << log_string; 
	}
*/
	template <class T>
	void		AddLog(const T& value)
	{
		_cache << value;
		if (_os!=NULL) (*_os) << value;
	}
	template <class T>
	void		AddLog(T& value)
	{
		_cache << value;
		if (_os!=NULL) (*_os) << value;
	}	
	// get full log string
	str_t		GetLog() { return _cache.str(); }
	// append the log to a file
	void 		SaveLog(bool clear=true)
	{
		if (_filename=="") throw "filename required in fjlib_log.h";
		std::ofstream of;
		of.open(_filename.c_str(),std::ios::out | std::ios::app);
		of << _cache.str();
		of.close();
		if (clear) _cache.str("");
	}
};	// end of class TFJLog

template <class T>
TFJLog& operator<<(TFJLog& log, const T& value)
{
	log.AddLog(value);
	return log;	
} 

template <class T>
TFJLog& operator<<(TFJLog& log, T& value)
{
	log.AddLog(value);
	return log;	
}

// shortcut for return symbol
#define cendl '\n'

}	// end of namespace

#endif
