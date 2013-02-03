#include "fjlib_memini.h"

namespace fjlib {

sect_find_pair TFJMemINI::find_section(const str_t& name)
{
	sectlist_type::iterator iter=_sects.find(name);
	if (iter!=_sects.end())
		return std::make_pair(true,sect_data_t(&(iter->second)));
	return std::make_pair(false,sect_data_t(NULL));
}

key_find_pair TFJMemINI::find_key(const str_t& key_name,
                                  const str_t& sect_name)
{
	sect_find_pair sf=find_section(sect_name);
	if (sf.first)
	{
		keylist_type::iterator iter=sf.second.keys().find(key_name);
		if (iter!=sf.second.keys().end())
			return std::make_pair(true,key_data_t(&(iter->second)));
	}
	return std::make_pair(false,key_data_t(NULL));
}

void TFJMemINI::trim_whitespaces(str_t& s)
{
	str_t szTrimChars=white_spaces;
	szTrimChars+=equal_indicator;
	int nPos, rPos;
	// trim left
	nPos = s.find_first_not_of(szTrimChars);
	if ( nPos > 0 ) s.erase(0, nPos);
	// trim right and return
	nPos = s.find_last_not_of(szTrimChars);
	rPos = s.find_last_of(szTrimChars);
	if ( rPos > nPos && rPos > -1)
    s.erase(rPos, s.size()-rPos);
}

str_t TFJMemINI::get_nextword(str_t& s)
{
	int nPos = s.find_first_of(equal_indicator);
	str_t sWord = str_t("");
	if ( nPos > -1 )
	{
		sWord = s.substr(0, nPos);
		s.erase(0, nPos+1);
	}
	else
	{
		sWord = s;
		s = str_t("");
	}
	trim_whitespaces(sWord);
	return sWord;
}

std::istream& operator>>(std::istream& in, TFJMemINI& c)
{
	c.clear();
	const size_t MAX_BUFFER_LEN=512;
	char buffer[MAX_BUFFER_LEN];
	str_t line,comment;
	str_t sect_name=str_t("");  // get global section
	bool_t done=false;
	while (!done)
	{
		memset(buffer, 0, MAX_BUFFER_LEN);
		in.getline(buffer, MAX_BUFFER_LEN);
		line=buffer;
		c.trim_whitespaces(line);
		done=(in.eof()||in.bad()||in.fail());
		if (line.size()==0) continue;
		if (line.find_first_of(c.comment_indicator)==0)
			comment+="\n"+line;
		else if (line.find_first_of('[')==0)
		{
			line.erase(0,1);
			line.erase(line.find_last_of(']'),1);
			sect_name=line;
			c.section(line).comment()=comment;
			comment=str_t("");
		}
		else
		{
			str_t key_name=c.get_nextword(line);
			str_t value=line;
			c.set_string(key_name,value,comment,sect_name);
			comment=str_t("");
		}
	}
	return in;
}

std::ostream& operator<<(std::ostream& out, const TFJMemINI& c)
{
	sectlist_type::const_iterator iter=c.begin();
	while (iter!=c.end())
	{
		out << "[" << iter->first << "]" << std::endl;
//		TFJINISectionData sect(&iter->second);
		const keylist_type* keys=&(iter->second.first);
		keylist_type::const_iterator it=keys->begin();
		while (it!=keys->end())
		{
			str_t comment=it->second.second;
			if (comment!="")
				out << "//" << comment << std::endl;
			out << it->first << "=" << it->second.first << std::endl;
			it++;
		}
		out << std::endl;
		iter++;
	}
	return out;
}


}  // end of namespace

