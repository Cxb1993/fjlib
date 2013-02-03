//---------------------------------------------------------------------------

#ifndef FJINIH
#define FJINIH
//---------------------------------------------------------------------------

#include "fjlib.h"
#include "fjlib_string.h"
#include "fjlib_map12.h"
//#include <sstream>

namespace fjlib {

//        typedef string_t str_t;
		typedef bool bool_t;
        typedef mapp<str_t,str_t,str_t> keylist_mapp;
        typedef keylist_mapp::map_type keylist_type;
        typedef mapp<str_t,keylist_type,str_t> sectlist_mapp;
        typedef sectlist_mapp::map_type sectlist_type;

struct TFJINISectionData {
        typedef sectlist_mapp::mapped_type mapped_type;
        typedef sectlist_mapp::data1_type data1_type;
        typedef sectlist_mapp::data2_type data2_type;
        mapped_type*    _sect;
        TFJINISectionData(mapped_type* v): _sect(v) {}
        TFJINISectionData& operator=(const TFJINISectionData& v)
        { _sect=v._sect; return *this; }
        data1_type&     keys() { return _sect->first; }
        data2_type&     comment() { return _sect->second; }
        data1_type&     keys() const { return _sect->first; }
        data2_type&     comment() const { return _sect->second; }
};

struct TFJINIKeyData {
        typedef keylist_mapp::mapped_type mapped_type;
        typedef keylist_mapp::data1_type data1_type;
        typedef keylist_mapp::data2_type data2_type;
        mapped_type*    _key;
        TFJINIKeyData(mapped_type* v): _key(v) {}
        TFJINIKeyData& operator=(const TFJINIKeyData& v)
        { _key=v._key; return *this; }
        data1_type&     value() { return _key->first; }
        data2_type&     comment() { return _key->second; }
        data1_type&     value() const { return _key->first; }
        data2_type&     comment() const { return _key->second; }
};

        const str_t CommentIndicator=str_t("//");
        const str_t EqualIndicator=str_t("=");
        const str_t WhiteSpaces =str_t(" \t\n\r");

        typedef TFJINISectionData sect_data_t;
        typedef TFJINIKeyData key_data_t;
		typedef std::pair<bool_t,sect_data_t> sect_find_pair;
		typedef std::pair<bool_t,key_data_t> key_find_pair;

class TFJMemINI {
private:
        sectlist_type   _sects;
        key_data_t      key(const str_t& name,
                            const sect_data_t& sect)
        { return key_data_t(&(sect.keys()[name])); }
public:
        TFJMemINI(): equal_indicator(EqualIndicator),
                     comment_indicator(CommentIndicator),
                     white_spaces(WhiteSpaces) {}
        ~TFJMemINI() {}

        str_t           equal_indicator;
        str_t           comment_indicator;
        str_t           white_spaces;

        // Section and key direct access

        sectlist_type::iterator
                        begin() { return _sects.begin(); }
        sectlist_type::iterator
                        end() { return _sects.end(); }
        sectlist_type::const_iterator
                        begin() const { return _sects.begin(); }
        sectlist_type::const_iterator
                        end() const { return _sects.end(); }

        sect_data_t     section(const str_t& name=str_t(""))
        { return sect_data_t(&(_sects[name])); }

        key_data_t      key(const str_t& key_name,
                            const str_t& sect_name=str_t(""))
        { return key(key_name,section(sect_name)); }

        // Section and key management

        sect_find_pair  find_section(const str_t& name=str_t(""));
        key_find_pair   find_key(const str_t& key_name,
                                 const str_t& sect_name=str_t(""));

        bool_t          del_section(const str_t& name=str_t(""))
        { return _sects.erase(name); }

        bool_t          del_key(const str_t& key_name,
                                const str_t& sect_name=str_t(""))
        {
			sect_find_pair sf=find_section(sect_name);
			if (sf.first) return sf.second.keys().erase(key_name);
			return false;
        }

        void            clear() { _sects.clear(); }
        size_t          sect_count() { return _sects.size(); }
        size_t          key_count(const str_t& sect_name)
        {
          sect_find_pair sf=find_section(sect_name);
          if (sf.first) return sf.second.keys().size();
          return 0;
        }

        // Data handling methods
        str_t           get_string(const str_t& _default,
                                   const str_t& key_name,
                                   const str_t& sect_name=str_t(""))
        {
          key_find_pair kf=find_key(key_name,sect_name);
          if (kf.first) return kf.second.value();
          return _default;
        }

		void			set_string(const str_t& key_name,
                                   const str_t& value,
                                   const str_t& comment=str_t(""),
                                   const str_t& sect_name=str_t(""))
        { key_data_t k=key(key_name,sect_name);
          k.value()=value; k.comment()=comment;
        }

        template <class T>
        T               get_stream(const T& _default,
                                   const str_t& key_name,
                                   const str_t& sect_name=str_t(""))
        {
			std::stringstream inout,tm;
			tm << _default;
			inout << get_string(tm.str(),key_name,sect_name).c_str();
			T value; inout >> value;
			return value;
        }

        template <class T>
        void            set_stream(const str_t& key_name,
                                   const T& value,
                                   const str_t& comment=str_t(""),
                                   const str_t& sect_name=str_t(""))
        {
			std::stringstream inout;
			inout << value;
			key_data_t k=key(key_name,sect_name);
			k.value()=inout.str();
			k.comment()=comment;
        }

        // Utility functions
        void            trim_whitespaces(str_t& s);
        str_t           get_nextword(str_t& s);
};

std::istream& operator>>(std::istream& in, TFJMemINI& c);
std::ostream& operator<<(std::ostream& out, const TFJMemINI& c);

}  // end of namespace


#endif
