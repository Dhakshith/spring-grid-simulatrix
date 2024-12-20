#include <memory>
#include <vector>
#include <map>
#include <string>
#include <iostream>

#define ISSPACE(x) ((x) == '\x20' || (x) == '\x09' || (x) == '\x0A' || (x) == '\x0D')

enum class JSONType {
	object,
	array,
	string,
	number
}; // I don't care about true, false, null, and decimals with a decimal point

class _JSON {
public:
	virtual JSONType type() const = 0;
	virtual ~_JSON(){};
	friend std::istream& operator>>(std::istream &in, _JSON json);
	friend std::ostream& operator<<(std::ostream &out, const _JSON& json);
};

class JSON {
public:
	std::unique_ptr<_JSON> data;
	JSON() {}
	friend std::istream& operator>>(std::istream &in, JSON& json);
	friend std::ostream& operator<<(std::ostream &out, const JSON& json);
};

class JSONNumber : public _JSON {
public:
	int n;
	bool valid;
	JSONNumber() : valid(false) {}
	JSONNumber(int n) : n(n), valid(true) {}
	JSONType type() const override {return JSONType::number;}
	friend std::istream& operator>>(std::istream &in, JSONNumber& jsonstr);
	friend std::ostream& operator<<(std::ostream &out, const JSONNumber& jsonn);
};

std::istream& operator>>(std::istream &in, JSONNumber& jsonn) {
	char c;
	jsonn.valid = false;
	while (in.read(&c, 1)) {
		if (!ISSPACE(c)) {
			if (!('0' <= c && c <= '9')) return in;
			else {
				jsonn.valid = true;
				break;
			}
		}
	}
	if (!jsonn.valid) return in;
	jsonn.n = c - '0';
	while (in.read(&c, 1)) {
		if ('0' <= c && c <= '9') jsonn.n = 10*jsonn.n + c-'0';
		else {
			in.putback(c);
			return in;
		}
	}
	jsonn.valid = false;
	return in;
}

std::ostream& operator<<(std::ostream &out, const JSONNumber& jsonn) {
	out << jsonn.n;
	return out;
}

class JSONString : public _JSON {
public:
	std::string str;
	bool valid;
	JSONString() : valid(false) {}
	JSONString(std::string str) : str(str), valid(true) {}
	JSONType type() const override {return JSONType::string;}
	friend std::istream& operator>>(std::istream &in, JSONString& jsonstr);
	friend std::ostream& operator<<(std::ostream &out, const JSONString& jsonstr);
};

std::istream& operator>>(std::istream &in, JSONString& jsonstr) {
	char c;
	jsonstr.valid = false;
	while (in.read(&c, 1)) {
		if (!ISSPACE(c)) {
			if (c != '"') return in;
			else {
				jsonstr.valid = true;
				break;
			}
		}
	}
	if (!jsonstr.valid) return in;
	while (in.read(&c, 1)) {
		if (c != '"') jsonstr.str += c; // I don't care about escaping and stuff
		else return in;
	}
	jsonstr.valid = false;
	return in;
}

std::ostream& operator<<(std::ostream &out, const JSONString& jsonstr) {
	out << '"' << jsonstr.str << '"';
	return out;
}

class JSONArray : public _JSON {
public:
	std::vector<JSON> objs;
	bool valid;
	JSONArray() : valid(false) {}
	JSONType type() const override {return JSONType::array;}
	void add(JSON&& json) {
		objs.push_back(std::forward<JSON>(json));
	}
	friend std::istream& operator>>(std::istream &in, JSONArray& jsonarr);
	friend std::ostream& operator<<(std::ostream &out, const JSONArray& jsonarr);
};

std::istream& operator>>(std::istream &in, JSONArray& jsonarr) {
	char c;
	jsonarr.valid = false;
	while (in.read(&c, 1)) {
		if (!ISSPACE(c)) {
			if (c != '[') return in;
			else {
				jsonarr.valid = true;
				break;
			}
		}
	}
	if (!jsonarr.valid) return in;
	
	jsonarr.valid = false;
	while(in.read(&c, 1)) {
		if (!ISSPACE(c)) {
			if (c == ']') {
				jsonarr.valid = true;
				return in;
			} else {
				jsonarr.valid = true;
				in.putback(c);
				break;
			}
		}
	}
	if (!jsonarr.valid) return in;

	while (true) {
		jsonarr.objs.emplace_back();
		in >> jsonarr.objs[jsonarr.objs.size()-1];
		if (!jsonarr.objs[jsonarr.objs.size()-1].data) {
			jsonarr.valid = false;
			jsonarr.objs.clear();
			return in;
		}
		while (in.read(&c, 1)) {
			if (!ISSPACE(c)) {
				if (c == ',') break;
				else if (c == ']') return in;
				else {
					jsonarr.valid = false;
					jsonarr.objs.clear();
					return in;
				}
			}
		}
	}
	jsonarr.valid = false;
	return in;
}

std::ostream& operator<<(std::ostream &out, const JSONArray& jsonarr) {
	out << '[';
	for (unsigned int i = 0; i < jsonarr.objs.size(); i++) {
		out << jsonarr.objs[i];
		if (i+1 != jsonarr.objs.size()) out << ", ";
	}
	out << ']';
	return out;
}

class JSONObj : public _JSON {
public:
	std::map<std::string, JSON> items;
	bool valid;
	JSONObj() {}
	JSONType type() const override {return JSONType::object;}
	void add(std::string key, JSON&& json) {
		items[key] = std::forward<JSON>(json);
	}
	friend std::istream& operator<<(std::istream &in, JSONObj& jsonobj);
	friend std::ostream& operator<<(std::ostream &out, const JSONObj& jsonobj);
};

std::istream& operator>>(std::istream &in, JSONObj& jsonobj) {
	char c;
	jsonobj.valid = false;
	while (in.read(&c, 1)) {
		if (!ISSPACE(c)) {
			if (c != '{') return in;
			else {
				jsonobj.valid = true;
				break;
			}
		}
	}
	if (!jsonobj.valid) return in;

	jsonobj.valid = false;
	while(in.read(&c, 1)) {
		if (!ISSPACE(c)) {
			if (c == '}') {
				jsonobj.valid = true;
				return in;
			} else if (c == '"') {
				jsonobj.valid = true;
				in.putback(c);
				break;
			} else return in;
		}
	}
	if (!jsonobj.valid) return in;

	while (true) {
		JSONString key;
		in >> key;
		if (!key.valid) {
			jsonobj.valid = false;
			jsonobj.items.clear();
			return in;
		}
		jsonobj.valid = false;
		while(in.read(&c, 1)) {
			if (!ISSPACE(c)) {
				if (c == ':') {
					jsonobj.valid = true;
					break;
				} else return in;
			}
		}
		if (!jsonobj.valid) {jsonobj.items.clear(); return in;}
		jsonobj.items.emplace(key.str, JSON());
		in >> jsonobj.items[key.str];
		if (!jsonobj.items[key.str].data) {
			jsonobj.valid = false;
			jsonobj.items.clear();
			return in;
		}
		jsonobj.valid = false;
		while (in.read(&c, 1)) {
			if (!ISSPACE(c)) {
				if (c == ',') break;
				else if (c == '}') {
					jsonobj.valid = true;
					return in;
				} else {jsonobj.items.clear(); return in;}
			}
		}
	}
	jsonobj.valid = false;
	jsonobj.items.clear();
	return in;
}

std::ostream& operator<<(std::ostream &out, const JSONObj& jsonobj) {
	out << '{';
	unsigned int x = 0;
	for (std::map<std::string, JSON>::const_iterator it = jsonobj.items.begin(); it != jsonobj.items.end(); it++) {
		out << '"' << it->first << '"' << ": " << it->second;
		if (x+1 != jsonobj.items.size()) out << ", ";
		x++;
	}
	out << '}';
	return out;
}

std::istream& operator>>(std::istream &in, JSON& json) {
	char c;
	while (in.read(&c, 1)) {
		if (ISSPACE(c)) continue;
		else if (c == '"') {
			in.putback(c);
			std::unique_ptr<JSONString> jsonstr(new JSONString);
			in >> *jsonstr;
			if (!jsonstr->valid) return in;
			json.data = std::move(jsonstr);
			return in;
		} else if (c == '[') {
			in.putback(c);
			std::unique_ptr<JSONArray> jsonarray(new JSONArray);
			in >> *jsonarray;
			if (!jsonarray->valid) return in;
			json.data = std::move(jsonarray);
			return in;
		} else if (c == '{') {
			in.putback(c);
			std::unique_ptr<JSONObj> jsonobj(new JSONObj);
			in >> *jsonobj;
			if (!jsonobj->valid) return in;
			json.data = std::move(jsonobj);
			return in;
		} else if ('0' <= c && c <= '9') {
			in.putback(c);
			std::unique_ptr<JSONNumber> jsonn(new JSONNumber);
			in >> *jsonn;
			if (!jsonn->valid) return in;
			json.data = std::move(jsonn);
			return in;
		} else {
			return in;
		}
	}
	return in;
}

std::ostream& operator<<(std::ostream &out, const JSON& json) {
	out << *(json.data);
	return out;
}

std::istream& operator>>(std::istream &in, _JSON& json) {
	if (json.type() == JSONType::object) in >> *((JSONObj*) &json);
	else if (json.type() == JSONType::array) in >> *((JSONArray*) &json);
	else if (json.type() == JSONType::string) in >> *((JSONString*) &json);
	else if (json.type() == JSONType::number) in >> *((JSONNumber*) &json);
	return in;
}

std::ostream& operator<<(std::ostream &out, const _JSON& json) {
	if (json.type() == JSONType::object) out << *((JSONObj*) &json);
	else if (json.type() == JSONType::array) out << *((JSONArray*) &json);
	else if (json.type() == JSONType::string) out << *((JSONString*) &json);
	else if (json.type() == JSONType::number) out << *((JSONNumber*) &json);
	return out;
}
