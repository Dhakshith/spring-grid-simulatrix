#include <fstream>
#include <cstdio>
#include "eigen3/Eigen/Core"
#include "json.cpp"

struct Stage {
	unsigned int n;
	std::vector<Eigen::Vector2d> home;
	std::vector<Eigen::Vector2d> fixedpoints;
	std::vector<std::vector<std::tuple<bool, int, double>>> connections;
};

/*
printf("%a\n", 16.0062511234);
	double tst = 0;
	std::string q;
	// std::cin >> q;
	char *end_;
	printf("%a\n", std::strtod(q.c_str(), &end_));
*/

bool getDouble(std::string s, double &d) {
	const char *cstr = s.c_str();
	char *t = NULL;

	d = std::strtod(cstr, &t);

	if (t != (cstr + s.size())) return false;

	return true;
}

std::string doubleToString(double d) {
	char s[32]; // I think 21 or so characters are actually needed
	sprintf(s, "%a", d);
	return std::string(s);
}

bool getVec2(JSONArray* vec, Eigen::Vector2d &v) {
	if (vec->objs.size() != 2) return false;
	if (vec->objs[0].data->type() != JSONType::string) return false;
	if (vec->objs[1].data->type() != JSONType::string) return false;
	
	JSONString *x = (JSONString*) vec->objs[0].data.get();
	JSONString *y = (JSONString*) vec->objs[1].data.get();

	if (!getDouble(x->str, v.x())) return false;
	if (!getDouble(y->str, v.y())) return false;

	return true;
}

void putStage(std::ostream &f, const Stage &stage) {

	std::unique_ptr<JSONNumber> jn(new JSONNumber(stage.n));
	std::unique_ptr<JSONArray> home(new JSONArray());
	std::unique_ptr<JSONArray> fixedpoints(new JSONArray());
	std::unique_ptr<JSONArray> connections(new JSONArray());

	for (unsigned int i = 0; i < stage.n; i++) {
		std::unique_ptr<JSONArray> vec(std::make_unique<JSONArray>());
		vec->objs.emplace_back();
		vec->objs.emplace_back();
		vec->objs[0].data = std::make_unique<JSONString>(doubleToString(stage.home[i].x()));
		vec->objs[1].data = std::make_unique<JSONString>(doubleToString(stage.home[i].y()));

		home->objs.emplace_back();
		home->objs[i].data = std::move(vec);
	}

	for (unsigned int i = 0; i < stage.fixedpoints.size(); i++) {
		std::unique_ptr<JSONArray> vec(std::make_unique<JSONArray>());
		vec->objs.emplace_back();
		vec->objs.emplace_back();
		vec->objs[0].data = std::make_unique<JSONString>(doubleToString(stage.fixedpoints[i].x()));
		vec->objs[1].data = std::make_unique<JSONString>(doubleToString(stage.fixedpoints[i].y()));

		fixedpoints->objs.emplace_back();
		fixedpoints->objs[i].data = std::move(vec);
	}

	for (unsigned int i = 0; i < stage.n; i++) {
		std::unique_ptr<JSONArray> conns(std::make_unique<JSONArray>());
		for (unsigned int j = 0; j < stage.connections[i].size(); j++) {
			std::unique_ptr<JSONArray> conn(std::make_unique<JSONArray>());
			conn->objs.emplace_back(); conn->objs.emplace_back(); conn->objs.emplace_back();
			conn->objs[0].data = std::make_unique<JSONNumber>(std::get<0>(stage.connections[i][j]));
			conn->objs[1].data = std::make_unique<JSONNumber>(std::get<1>(stage.connections[i][j]));
			conn->objs[2].data = std::make_unique<JSONString>(doubleToString(std::get<2>(stage.connections[i][j])));
			conns->objs.emplace_back();
			conns->objs[j].data = std::move(conn);
		}

		connections->objs.emplace_back();
		connections->objs[i].data = std::move(conns);
	}


	std::unique_ptr<JSONObj> jobj(new JSONObj());
	jobj->items.emplace("n", JSON()); jobj->items["n"].data = std::move(jn);
	jobj->items.emplace("home", JSON()); jobj->items["home"].data = std::move(home);
	jobj->items.emplace("fixedpoints", JSON()); jobj->items["fixedpoints"].data = std::move(fixedpoints);
	jobj->items.emplace("connections", JSON()); jobj->items["connections"].data = std::move(connections);
	
	JSON j; j.data = std::move(jobj);

	f << j;
}

void putStage(std::string fname, const Stage &stage) {
	std::ofstream f(fname);
	putStage(f, stage);
}

bool getStage(std::istream &f, Stage &stage) {
	JSON j;

	f >> j;

	if (!j.data) return false;
	
	char c;
	while (f.read(&c, 1)) {
		if (!ISSPACE(c)) return false;
	}

	if (j.data->type() != JSONType::object) return false;
	
	JSONObj *jobj = (JSONObj*) j.data.get();

	if (jobj->items.find("n") == jobj->items.end()) return false;
	if (jobj->items.find("home") == jobj->items.end()) return false;
	if (jobj->items.find("fixedpoints") == jobj->items.end()) return false;
	if (jobj->items.find("connections") == jobj->items.end()) return false;

	if (jobj->items["n"].data->type() != JSONType::number) return false;
	if (jobj->items["home"].data->type() != JSONType::array) return false;
	if (jobj->items["fixedpoints"].data->type() != JSONType::array) return false;
	if (jobj->items["connections"].data->type() != JSONType::array) return false;

	JSONNumber *n_ = (JSONNumber*) jobj->items["n"].data.get();
	JSONArray *home_ = (JSONArray*) jobj->items["home"].data.get();
	JSONArray *fixedpoints_ = (JSONArray*) jobj->items["fixedpoints"].data.get();
	JSONArray *connections_ = (JSONArray*) jobj->items["connections"].data.get();

	if (n_->n < 0) return false;
	stage.n = n_->n;
	
	if (home_->objs.size() != stage.n) return false;
	if (connections_->objs.size() != stage.n) return false;

	stage.home.clear(); stage.fixedpoints.clear(); stage.connections.clear();
	stage.home.assign(stage.n, Eigen::Vector2d(0, 0)); stage.connections.assign(stage.n, std::vector<std::tuple<bool, int, double>>());
	stage.fixedpoints.assign(fixedpoints_->objs.size(), Eigen::Vector2d(0, 0));

	for (unsigned int i = 0; i < stage.n; i++) {
		if (home_->objs[i].data->type() != JSONType::array) return false;
		if (!getVec2((JSONArray*) home_->objs[i].data.get(), stage.home[i])) return false;
	}

	for (unsigned int i = 0; i < fixedpoints_->objs.size(); i++) {
		if (fixedpoints_->objs[i].data->type() != JSONType::array) return false;
		if (!getVec2((JSONArray*) fixedpoints_->objs[i].data.get(), stage.fixedpoints[i])) return false;
	}

	for (unsigned int i = 0; i < stage.n; i++) {
		if (connections_->objs[i].data->type() != JSONType::array) return false;
		JSONArray *el__ = (JSONArray*) connections_->objs[i].data.get();
		for (unsigned int j = 0; j < el__->objs.size(); j++) {
			if (el__->objs[j].data->type() != JSONType::array) return false;
			JSONArray *conn = (JSONArray*) el__->objs[j].data.get();
			if (conn->objs.size() != 3) return false;
			if (conn->objs[0].data->type() != JSONType::number) return false;
			if (conn->objs[1].data->type() != JSONType::number) return false;
			if (conn->objs[2].data->type() != JSONType::string) return false;
			int a = ((JSONNumber*)conn->objs[0].data.get())->n;
			int b = ((JSONNumber*)conn->objs[1].data.get())->n;
			std::string s = ((JSONString*)conn->objs[2].data.get())->str;
			if (a != 0 && a != 1) return false;
			if (a == 0) {
				if (!(0 <= b && (unsigned int) b < stage.n)) return false;
			} else {
				if (!(0 <= b && (unsigned int) b < stage.fixedpoints.size())) return false;
			}
			double d;
			if (!getDouble(s, d)) return false;
			if (d < 0) return false;
			stage.connections[i].emplace_back(a, b, d);
		}
	}

	return true;
}

bool getStage(std::string fname, Stage &stage) {
	std::ifstream f(fname);
	return getStage(f, stage);
}

/*
int main() {
	Stage test;
	std::cout << getStage("Sample.json", test) << '\n';
	std::cout << test.n << '\n';
	for (auto u : test.home) std::cout << u.x() << ':' << u.y() << ' '; std::cout << '\n';
	for (auto u : test.fixedpoints) std::cout << u.x() << ':' << u.y() << ' '; std::cout << '\n';
	for (auto v : test.connections) {
		for (auto u : v) std::cout << std::get<0>(u) << ' ' << std::get<1>(u) << ' ' << std::get<2>(u) << ';';
		std::cout << '\n';
	} std::cout << '\n';

	putStage("SampleOp.json", test);
	std::cout << getStage("SampleOp.json", test) << '\n';

	for (auto u : test.home) std::cout << u.x() << ':' << u.y() << ' '; std::cout << '\n';
	for (auto u : test.fixedpoints) std::cout << u.x() << ':' << u.y() << ' '; std::cout << '\n';
	for (auto v : test.connections) {
		for (auto u : v) std::cout << std::get<0>(u) << ' ' << std::get<1>(u) << ' ' << std::get<2>(u) << ';';
		std::cout << '\n';
	} std::cout << '\n';

	return 0;
}
*/
