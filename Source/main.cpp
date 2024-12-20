#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#include <emscripten/html5.h>
#endif

#include <iostream>
#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <functional>
#include <algorithm>
#include <cmath>

#include "eigen3/Eigen/Core"
#include "shapes.cpp"
#include "widgets.cpp"
#include "solver.cpp"
#include "stage.cpp"

const double PI = 3.141592653589793238463;

// I don't want to do portrait. Let it be all landscape.
const int scr_width_landscape = 800;
const int scr_height_landscape = 600;
const int scr_width_portrait = 600;
const int scr_height_portrait = 800;
const double aspect_ratio_landscape = (double) scr_width_landscape / scr_height_landscape;
const double aspect_ratio_portrait = (double) scr_width_portrait / scr_height_portrait;

std::string globaljson;
bool received;

class Scene {
public:
	virtual void init() = 0;
	virtual void place() = 0;
    virtual void eventhandler(SDL_Event& e) = 0;
    virtual void run(double delta) = 0;

    virtual ~Scene(){};
};

struct Globals {
	int screen_width;
	int screen_height;
	bool portrait;
	SDL_Rect rect;
	double scaleRatio;

	SDL_Window *window;
	SDL_Renderer *renderer;
	TTF_Font *font1, *font2;
	SDL_RendererInfo rendInfo;

	std::vector<Label> *labels;
	std::vector<Button> *buttons;
	std::vector<Slider> *sliders;
	std::vector<Numberinput> *numberinputs;
	std::vector<Checkbox> *chkbxes;

	Scene *scene;
	bool quit;
	char currchar;
	Uint64 last_milli;
	Numberinput* inputbox;

	void changescene(Scene*);
};

int translatedX(Globals *globals, int x) {
	return globals->rect.x + x * globals->scaleRatio;
}

int translatedY(Globals *globals, int y) {
	return globals->rect.y + y * globals->scaleRatio;
}

std::pair<int, int> translatedCoords(Globals *globals, int x, int y) {
	return std::pair<int, int>(globals->rect.x + x * globals->scaleRatio, globals->rect.y + y * globals->scaleRatio);
}

SDL_Point reverseTranslate(Globals *globals, int x, int y) {
	return (SDL_Point) {(int) ((x - globals->rect.x) / globals->scaleRatio), (int) ((y - globals->rect.y) / globals->scaleRatio)};
}

void Globals::changescene(Scene *newscene) {
	labels->clear();
	buttons->clear();
	sliders->clear();
	numberinputs->clear();
	chkbxes->clear();

	delete scene;
	scene = newscene;
	scene->init();
	scene->place();
}

void fill_screen(void *userData) {
	double w, h;
	emscripten_get_element_css_size("canvas", &w, &h);

	Globals* globals = (Globals*) userData;
	globals->screen_width = (int) w; globals->screen_height = (int) h;

	//globals->portrait = globals->screen_width < globals->screen_height;
	globals->portrait = false;

	{
		double aspect_ratio = globals->portrait ? aspect_ratio_portrait : aspect_ratio_landscape;
		if (globals->screen_width >= globals->screen_height * aspect_ratio) {
			globals->rect.h = globals->screen_height;
			globals->rect.w = globals->rect.h * aspect_ratio;
		} else {
			globals->rect.w = globals->screen_width;
			globals->rect.h = globals->rect.w / aspect_ratio;
		}

		globals->rect.x = (globals->screen_width - globals->rect.w)/2;
		globals->rect.y = (globals->screen_height - globals->rect.h)/2;
	}

	globals->scaleRatio = (double) globals->rect.w / (globals->portrait ? scr_width_portrait : scr_width_landscape);

	TTF_SetFontSize(globals->font1, 15 * globals->scaleRatio);
	TTF_SetFontSize(globals->font2, 24 * globals->scaleRatio);
	if (globals->scene) globals->scene->place();

	SDL_SetWindowSize((SDL_Window*) globals->window, globals->screen_width, globals->screen_height);
}

EM_BOOL window_resized_callback(int eventType, const void *reserved, void *userData) {
	fill_screen(userData);
	
	return true;
}

EM_BOOL orientation_change(int eventType, const EmscriptenOrientationChangeEvent *orientationChangeEvent, void *userData) {
	fill_screen(userData);

	return true;
}

bool initialize(struct Globals *globals) {
	if (SDL_Init(SDL_INIT_VIDEO) < 0) {
		printf("Could not initialize SDL: %s\n", SDL_GetError());
		return true;
	}

	if (TTF_Init() < 0) {
	   printf("TTF Init Error: %s\n", TTF_GetError());
	   return true;
	}

	globals->window = SDL_CreateWindow("Spring Grid Simulator", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, 800, 600, SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);
	if (!globals->window) {
		printf("Could not create window: %s\n", SDL_GetError());
		return true;
	}

	globals->renderer = SDL_CreateRenderer(globals->window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
	if (!globals->renderer) {
		printf("Could not create renderer: %s\n", SDL_GetError());
		return true;
	}

	globals->font1 = TTF_OpenFont("Resources/Inconsolata-Thin.ttf", 15);
	if (!globals->font1) {
		printf("Font can't be opened: %s\n", SDL_GetError());
		return true;
	}

	globals->font2 = TTF_OpenFont("Resources/Inconsolata-Thin.ttf", 24);
	if (!globals->font2) {
		printf("Font can't be opened: %s\n", SDL_GetError());
		return true;
	}

	if (SDL_GetRendererInfo(globals->renderer, &globals->rendInfo) < 0) {
		printf("Cannot get renderer info: %s\n", SDL_GetError());
		return true;
	}

	// Soft Fullscreen
	EmscriptenFullscreenStrategy strat;
	strat.scaleMode = EMSCRIPTEN_FULLSCREEN_SCALE_STRETCH;
	strat.filteringMode = EMSCRIPTEN_FULLSCREEN_FILTERING_DEFAULT;
	strat.canvasResizedCallback = window_resized_callback;
	strat.canvasResizedCallbackUserData = globals;
	emscripten_enter_soft_fullscreen("canvas", &strat);

	emscripten_set_orientationchange_callback(globals, true, orientation_change);

	return false;
}

class D2Logic {
public:
	int n;
	Eigen::MatrixXd xycoeff;
	Eigen::VectorXd constfacx;
	Eigen::VectorXd constfacy;

	D2Logic() {}
	D2Logic(Stage stage) : n(stage.n) {
		if (n != 0) {
			xycoeff = Eigen::MatrixXd::Zero(n, n);
			constfacx = Eigen::VectorXd::Zero(n);
			constfacy = Eigen::VectorXd::Zero(n);
			for (int i = 0; i < n; i++) {
				for (auto j : stage.connections[i]) {
					if (std::get<0>(j) == 0) {
						xycoeff(i, std::get<1>(j)) = std::get<2>(j);
					} else {
						constfacx(i) += stage.fixedpoints[std::get<1>(j)].x() * std::get<2>(j);
						constfacy(i) += stage.fixedpoints[std::get<1>(j)].y() * std::get<2>(j);
					}

					xycoeff(i, i) -= std::get<2>(j);
				}
			}
		}
	}
};

class Menu : public Scene {
public:
	struct Globals *globals;

	Menu(struct Globals*);

	void init() override;
	void place() override;
	void eventhandler(SDL_Event&) override;
	void run(double) override;
};

enum class SelectType {
	None,
	Movingpt,
	Fixedpt,
	Conn, Conn2,
	LineM, LineM2,
	LineF, LineF2,
	CircleM, CircleM2,
	CircleF, CircleF2,
	Remove
};

class Create : public Scene {
public:
	struct Globals *globals;

	Stage stage;
	inline static int RAD = 10;
	inline static int fRAD = 5;
	int sqside;

	SelectType selected;
	struct {
		bool underType;
		int under;
	} selectedPrev;

	bool beingdraggedFixed, beingdraggedOrNot;
	int beingdragged;
	bool highlighted;
	int highlight;

	int linmassinp, linfixinp, circlemassinp, circlefixinp;

	Create(struct Globals*);
	Create(struct Globals*, Stage);

	void init() override;
	void place() override;
	void removefixedpoint(unsigned int);
	void removemovingpoint(unsigned int);
	void eventhandler(SDL_Event&) override;
	void run(double) override;
};

void placeButtonAndTakeCareOfEverything(Button &b, Globals* globals, int bwidth, int x, int y, bool useCoordAsCenter=false) {
	b = Button(b.renderer, b.text, b.onclick, bwidth*globals->scaleRatio, translatedCoords(globals, x, y), b.font, useCoordAsCenter, b.depressed, b.enabled);
}

class Mainscene : public Scene {
public:
	struct Globals *globals;

	Stage stage;
	int n;
	inline static int RAD = 10;
	inline static int velRAD = 4;
	inline static double DAMPMAX = 0.5;
	double DAMP;
	int sqside;

	std::vector<Eigen::Vector2d> pos;
	std::vector<Eigen::Vector2d> vel;
	std::vector<std::vector<bool>> collnstate;
	D2Logic d2logic;
	Solver xsolver, ysolver;

	bool PAUSED, SHOWVELOCITY, COLLIDE, changeoccurred;
	int pausebtn, resetposbtn, resetvelbtn, menubtn;
	int dampening;
	int showvel, collisions;

	bool beingdraggedVelocity, beingdraggedOrNot;
	int beingdragged;

	double currtime;

	Mainscene(struct Globals *globals, Stage stage) : globals(globals), stage(stage), n(stage.n), DAMP(0.0), sqside(600), pos(stage.home), vel(n, Eigen::Vector2d(0, 0)), collnstate(n, std::vector<bool>(n, false)) {
		if (n != 0) {
			d2logic = D2Logic(stage);
			xsolver = Solver(n, d2logic.xycoeff, DAMP, d2logic.constfacx);
			ysolver = Solver(n, d2logic.xycoeff, DAMP, d2logic.constfacy);
		}
	}

	void solve() {
		if (n == 0) return;
		Eigen::VectorXd x0(n), y0(n), vx0(n), vy0(n);

		for (int i = 0; i < n; i++) {
			x0[i] = pos[i].x();
			y0[i] = pos[i].y();
			vx0[i] = vel[i].x();
			vy0[i] = vel[i].y();
		}

		xsolver.solve(x0, vx0, 0);
		ysolver.solve(y0, vy0, 0);
	}

	void init() override {
		PAUSED = true;
		SHOWVELOCITY = false;
		COLLIDE = false;
		changeoccurred = false;
		beingdraggedVelocity = false;
		beingdraggedOrNot = false;
		currtime = 0;

		if (n != 0) solve();

		globals->buttons->emplace_back(globals->renderer, "Resume", [this] (Button &btn) {
			btn.updateText(btn.text == "Pause" ? "Resume" : "Pause");
			if (this->PAUSED && this->changeoccurred) {
				this->changeoccurred = false;
				this->solve();
				this->currtime = 0;
			}

			this->PAUSED = !this->PAUSED;
			(*this->globals->buttons)[this->resetposbtn].setenabled(this->PAUSED);
			(*this->globals->buttons)[this->resetvelbtn].setenabled(this->PAUSED);
			(*this->globals->sliders)[this->dampening].setenabled(this->PAUSED);
			(*this->globals->chkbxes)[this->collisions].setenabled(this->PAUSED);
		}, 80, std::make_pair<int, int>(sqside+5, 20), globals->font1, false);
		pausebtn = globals->buttons->size() - 1;
		
		this->globals->buttons->emplace_back(globals->renderer, "Reset position", [this] (Button &b) {
			this->pos = this->stage.home;
			this->changeoccurred = true;
		}, 80, std::make_pair<int, int>(sqside+5, 45), globals->font1, false);
		resetposbtn = globals->buttons->size() - 1;

		globals->buttons->emplace_back(globals->renderer, "Zero velocity", [this] (Button &b) {
			this->vel = std::vector<Eigen::Vector2d>(n, Eigen::Vector2d(0, 0));
			this->changeoccurred = true;
		}, 80, std::make_pair<int, int>(sqside+5, 70), globals->font1, false);
		resetvelbtn = globals->buttons->size() - 1;

		globals->buttons->emplace_back(globals->renderer, "Menu", [=] (Button &b) {globals->changescene(new Menu(globals));}, 80, std::make_pair<int, int>(sqside+5, 190), globals->font1, false);
		menubtn = globals->buttons->size() - 1;
		
		globals->buttons->emplace_back(globals->renderer, "Modify", [this] (Button &b) {globals->changescene(new Create(globals, this->stage));}, 80, std::make_pair<int, int>(sqside+5, 220), globals->font1, false);

		globals->sliders->emplace_back(globals->renderer, 80, 10, std::make_pair<int, int>(sqside+5, 95), [this] (Slider &sld) {
			this->DAMP = this->xsolver.damp = this->ysolver.damp = this->DAMPMAX * sld.value;
			this->changeoccurred = true;
		});
		dampening = globals->sliders->size() - 1;

		globals->chkbxes->emplace_back(globals->renderer, "Show velocity", [this] (Checkbox &chkbx) {
			this->SHOWVELOCITY = chkbx.ticked;
		}, std::pair<int, int>(sqside+5, 120), globals->font1);
		showvel = globals->chkbxes->size() - 1;

		globals->chkbxes->emplace_back(globals->renderer, "Collisions", [this] (Checkbox &chkbx) {
			this->COLLIDE = chkbx.ticked;
		}, std::pair<int, int>(sqside+5, 145), globals->font1);
		collisions = globals->chkbxes->size() - 1;
		
		(*globals->sliders)[dampening].value = DAMP / DAMPMAX;

		beingdraggedVelocity = false;
		beingdraggedOrNot = false;
	}

	void place() override {
		if (globals->portrait) {
		} else {
			placeButtonAndTakeCareOfEverything((*globals->buttons)[0], globals, 80, sqside+5, 20);
			placeButtonAndTakeCareOfEverything((*globals->buttons)[1], globals, 80, sqside+5, 45);
			placeButtonAndTakeCareOfEverything((*globals->buttons)[2], globals, 80, sqside+5, 70);
			placeButtonAndTakeCareOfEverything((*globals->buttons)[3], globals, 80, sqside+5, 190);
			placeButtonAndTakeCareOfEverything((*globals->buttons)[4], globals, 80, sqside+5, 220);
			(*globals->chkbxes)[0] = Checkbox(globals->renderer, (*globals->chkbxes)[0].text, (*globals->chkbxes)[0].onchange, translatedCoords(globals, sqside+5, 120), globals->font1, (*globals->chkbxes)[0].ticked, (*globals->chkbxes)[0].depressed, (*globals->chkbxes)[0].enabled);
			(*globals->chkbxes)[1] = Checkbox(globals->renderer, (*globals->chkbxes)[1].text, (*globals->chkbxes)[1].onchange, translatedCoords(globals, sqside+5, 145), globals->font1, (*globals->chkbxes)[1].ticked, (*globals->chkbxes)[1].depressed, (*globals->chkbxes)[1].enabled);
			(*globals->sliders)[0] = Slider(globals->renderer, 80*globals->scaleRatio, 10*globals->scaleRatio, translatedCoords(globals, sqside+5, 95), (*globals->sliders)[0].onchange, (*globals->sliders)[0].colorleft, (*globals->sliders)[0].colorright, (*globals->sliders)[0].colorptr, (*globals->sliders)[0].depressed, (*globals->sliders)[0].enabled);
		}
	}

	void updatemasses(double t) {
		if (n == 0) return;
		auto [xt, vxt] = xsolver.soln(t);
		auto [yt, vyt] = ysolver.soln(t);
		
		for (int i = 0; i < n; i++) {
			pos[i].x() = xt[i];
			pos[i].y() = yt[i];
			vel[i].x() = vxt[i];
			vel[i].y() = vyt[i];
		}
	}

	void eventhandler(SDL_Event& e) override {
		int x, y;
		SDL_GetMouseState(&x, &y);
		auto ourCoordXY = reverseTranslate(globals, x, y); x = ourCoordXY.x; y = ourCoordXY.y;
		Eigen::Vector2d mousepos(x, y);
		
		if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT) {
			if (PAUSED == true) {
				for (int i = 0; i < n; i++) {
					if ((mousepos-sqside*pos[i]).norm()<=RAD) {
						beingdragged = i;
						beingdraggedVelocity = false;
						beingdraggedOrNot = true;
						changeoccurred = true;
					}

					if (SHOWVELOCITY && (((pos[i]+vel[i])*sqside)-mousepos).norm() <= velRAD*globals->scaleRatio) {
						beingdragged = i;
						beingdraggedVelocity = true;
						beingdraggedOrNot = true;
						changeoccurred = true;
					}
				}
			}
		}

		if (e.type == SDL_MOUSEBUTTONUP && e.button.button == SDL_BUTTON_LEFT) {
			beingdraggedOrNot = false;
		}

		if (e.type == SDL_MOUSEMOTION) {
			if (beingdraggedOrNot) {
				if (beingdraggedVelocity) {
					vel[beingdragged] = (mousepos/sqside)-pos[beingdragged];
				} else {
					pos[beingdragged] = mousepos/sqside;
				}
			}
		}

		if (e.type == SDL_KEYDOWN && ' ' == e.key.keysym.sym)
			(*globals->buttons)[pausebtn].onclick((*globals->buttons)[pausebtn]);
	}

	void run(double delta) override {
		if (n == 0) return;
		if (!PAUSED) {
			updatemasses(currtime);
			currtime += 0.495 * delta/36;

			if (COLLIDE) {
				bool colln = false;
				for (int m1 = 0; m1 < n; m1++) {
					for (int m2 = m1+1; m2 < n; m2++) {
						if (sqside * (pos[m1]-pos[m2]).norm() <= RAD*2) {
							if (collnstate[m1][m2]) continue;
							collnstate[m1][m2] = true;
							Eigen::Vector2d ncap = pos[m1] - pos[m2];
							if (ncap.norm() != 0) {
								ncap = ncap / ncap.norm();
								if ((vel[m1] - vel[m2]).dot(ncap) < 0) {
									Eigen::Vector2d finalvelm1 = vel[m1] + ncap * (ncap.dot(vel[m2] - vel[m1]));
									vel[m2] = vel[m2] + ncap * (ncap.dot(vel[m1] - vel[m2]));
									vel[m1] = finalvelm1;
									colln = true;
								}
							}
						} else {
							collnstate[m1][m2] = false;
						}
					}
				}

				if (colln) {solve(); currtime = 0;}
			}
		}

		for (auto p : pos) {
			auto [tmpX, tmpY] = translatedCoords(globals, sqside*p.x(), sqside*p.y());
			Circle(globals->renderer, tmpX, tmpY, RAD*globals->scaleRatio, (SDL_Color) {0, 0, 0}, true);
		}

		if (SHOWVELOCITY) {
			for (int i = 0; i < n; i++) {
				Eigen::Vector2d _1_ = pos[i], _2_ = vel[i];
				SDL_SetRenderDrawColor(globals->renderer, 0x00, 0x00, 0xff, 0xff);
				auto [tmpX1, tmpY1] = translatedCoords(globals, sqside*_1_.x(), sqside*_1_.y());
				auto [tmpX2, tmpY2] = translatedCoords(globals, sqside*(_1_+_2_).x(), sqside*(_1_+_2_).y());
				SDL_RenderDrawLineF(globals->renderer, tmpX1, tmpY1, tmpX2, tmpY2);
				Circle(globals->renderer, tmpX2, tmpY2, velRAD*globals->scaleRatio, (SDL_Color) {0, 0, 255}, true);
			}
		}
		
		for (int c = 0; c < n; c++) {
			for (auto [typ, endpt_, weight] : stage.connections[c]) {
				Eigen::Vector2d startpt = pos[c];
				Eigen::Vector2d endpt = typ == 0 ? pos[endpt_] : stage.fixedpoints[endpt_];
				SDL_SetRenderDrawColor(globals->renderer, 0x00, 0x00, 0x00, 0xff);
				auto [tmpX1, tmpY1] = translatedCoords(globals, sqside*startpt.x(), sqside*startpt.y());
				auto [tmpX2, tmpY2] = translatedCoords(globals, sqside*endpt.x(), sqside*endpt.y());
				SDL_RenderDrawLineF(globals->renderer, tmpX1, tmpY1, tmpX2, tmpY2);
			}
		}
	}
};

Stage squarestage(unsigned int n, double strength=0.2) {
	Stage stage{n*n, std::vector<Eigen::Vector2d>(n*n, Eigen::Vector2d(0, 0)), std::vector<Eigen::Vector2d>((n+2)*(n+2), Eigen::Vector2d(0, 0)), std::vector<std::vector<std::tuple<bool, int, double>>>(n*n, std::vector<std::tuple<bool, int, double>>())};

	for (unsigned int i = 0; i < n; i++) {
		for (unsigned int j = 0; j < n; j++) {
			stage.home[i*n+j].x() = (double)(i+1)/(n+1);
			stage.home[i*n+j].y() = (double)(j+1)/(n+1);
		}
	}

	for (unsigned int i = 0; i < n+2; i++) {
		for (unsigned int j = 0; j < n+2; j++) {
			stage.fixedpoints[i*(n+2)+j].x() = (double)i/(n+1);
			stage.fixedpoints[i*(n+2)+j].y() = (double)j/(n+1);
		}
	}

	for (unsigned int i = 0; i < n; i++) {
		for (unsigned int j = 0; j < n; j++) {
			constexpr std::array<std::pair<int, int>, 4> dirs{std::pair<int, int>(0, -1), std::pair<int, int>(0, 1), std::pair<int, int>(-1, 0), std::pair<int, int>(1, 0)};
			for (auto [u, v] : dirs) {
				stage.connections[j*n+i].push_back((0 <= i+u && i+u < n && 0 <= j+v && j+v < n) ? std::tuple<bool, int, double>(false, (j+v)*n+(i+u), strength) : std::tuple<bool, int, double>(true, (j+v+1)*(n+2)+(i+u+1), strength));
			}
		}
	}

	putStage("s2.json", stage);

	return stage;
}

class Enter : public Scene {
public:
	struct Globals *globals;

	Enter(struct Globals *globals) : globals(globals) {}

	void init() override {
		globals->labels->emplace_back(globals->renderer, "Select size: ", std::make_pair<int, int>(50, 50), globals->font2, (SDL_Color) {0, 0, 0}, false);

		globals->buttons->emplace_back(globals->renderer, "Open file", [=] (Button &b) {
			Stage s{};
			if (getStage("s2.json", s))
				globals->changescene(new Mainscene(globals, s));
		}, 10, std::make_pair<int, int>(50, 100), globals->font1, false);
	
		for (int i = 1; i <= 10; i++)
			globals->buttons->emplace_back(globals->renderer, std::to_string(i), [=] (Button &b) {globals->changescene(new Mainscene(globals, squarestage(i)));}, 10, std::make_pair<int, int>(180 + i*20, 52), globals->font1, false);
	}

	void place() override {
		if (globals->portrait) {
			(*globals->labels)[0] = Label(globals->renderer, (*globals->labels)[0].text, translatedCoords(globals, 50, 50), globals->font2, (SDL_Color) {0, 0, 0}, false);
			placeButtonAndTakeCareOfEverything((*globals->buttons)[0], globals, 10, 50, 100);
			for (int i = 1; i <= 10; i++)
				placeButtonAndTakeCareOfEverything((*globals->buttons)[i], globals, 10, 240+i*20, 58);
		} else {
			(*globals->labels)[0] = Label(globals->renderer, (*globals->labels)[0].text, translatedCoords(globals, 50, 50), globals->font2, (SDL_Color) {0, 0, 0}, false);
			placeButtonAndTakeCareOfEverything((*globals->buttons)[0], globals, 10, 50, 100);
			for (int i = 1; i <= 10; i++)
				placeButtonAndTakeCareOfEverything((*globals->buttons)[i], globals, 10, 240+i*20, 47);
		}
	}

	void eventhandler(SDL_Event& e) override {}
	void run(double delta) override {}
};

Create::Create(struct Globals *globals) : globals(globals), stage(), sqside(globals->screen_height), selected(SelectType::None), beingdraggedFixed(true), beingdraggedOrNot(false), highlighted(false) {}

Create::Create(struct Globals *globals, Stage stage) : globals(globals), stage(stage), sqside(globals->screen_height), selected(SelectType::None), beingdraggedFixed(true), beingdraggedOrNot(false), highlighted(false) {}

extern "C" {
	void setJsonString(char *j) {
		globaljson = std::string(j);
		free(j);
		received = true;
	}
}

void Create::init() {
	globals->buttons->emplace_back(globals->renderer, "Menu", [=] (Button &btn) {globals->changescene(new Menu(globals));}, 80, std::make_pair<int, int>(sqside+5, 20), globals->font1, false);

	globals->buttons->emplace_back(globals->renderer, "Add mass", [this] (Button &btn) {this->selected = SelectType::Movingpt;}, 80, std::make_pair<int, int>(sqside+5, 45), globals->font1, false);

	globals->buttons->emplace_back(globals->renderer, "Add fixed point", [this] (Button &btn) {this->selected = SelectType::Fixedpt;}, 80, std::make_pair<int, int>(sqside+5, 70), globals->font1, false);

	globals->buttons->emplace_back(globals->renderer, "Add connection", [this] (Button &btn) {this->selected = SelectType::Conn;}, 80, std::make_pair<int, int>(sqside+5, 95), globals->font1, false);

	globals->buttons->emplace_back(globals->renderer, "Line mass", [this] (Button &btn) {this->selected = SelectType::LineM;}, 110, std::make_pair<int, int>(sqside+5, 120), globals->font1, false);
	globals->numberinputs->emplace_back(globals->renderer, "", 30, std::make_pair<int, int>(sqside+120, 120), globals->font1, 2);
	linmassinp = globals->numberinputs->size() - 1;
	globals->buttons->emplace_back(globals->renderer, "\u2191", [=] (Button &btn) {(*globals->numberinputs)[linmassinp].inc();}, 10, std::make_pair<int, int>(sqside+155, 120), globals->font1, false);
	globals->buttons->emplace_back(globals->renderer, "\u2193", [=] (Button &btn) {(*globals->numberinputs)[linmassinp].dec();}, 10, std::make_pair<int, int>(sqside+170, 120), globals->font1, false);

	globals->buttons->emplace_back(globals->renderer, "Line points", [this] (Button &btn) {this->selected = SelectType::LineF;}, 110, std::make_pair<int, int>(sqside+5, 145), globals->font1, false);
	globals->numberinputs->emplace_back(globals->renderer, "", 30, std::make_pair<int, int>(sqside+120, 145), globals->font1, 2);
	linfixinp = globals->numberinputs->size() - 1;
	globals->buttons->emplace_back(globals->renderer, "\u2191", [=] (Button &btn) {(*globals->numberinputs)[linfixinp].inc();}, 10, std::make_pair<int, int>(sqside+155, 145), globals->font1, false);
	globals->buttons->emplace_back(globals->renderer, "\u2193", [=] (Button &btn) {(*globals->numberinputs)[linfixinp].dec();}, 10, std::make_pair<int, int>(sqside+170, 145), globals->font1, false);

	globals->buttons->emplace_back(globals->renderer, "Circle mass", [this] (Button &btn) {this->selected = SelectType::CircleM;}, 110, std::make_pair<int, int>(sqside+5, 170), globals->font1, false);
	globals->numberinputs->emplace_back(globals->renderer, "", 30, std::make_pair<int, int>(sqside+120, 170), globals->font1, 2);
	circlemassinp = globals->numberinputs->size() - 1;
	globals->buttons->emplace_back(globals->renderer, "\u2191", [this] (Button &btn) {(*globals->numberinputs)[circlemassinp].inc();}, 10, std::make_pair<int, int>(sqside+155, 170), globals->font1, false);
	globals->buttons->emplace_back(globals->renderer, "\u2193", [this] (Button &btn) {(*globals->numberinputs)[circlemassinp].dec();}, 10, std::make_pair<int, int>(sqside+170, 170), globals->font1, false);

	globals->buttons->emplace_back(globals->renderer, "Circle points", [this] (Button &btn) {this->selected = SelectType::CircleF;}, 110, std::make_pair<int, int>(sqside+5, 195), globals->font1, false);
	globals->numberinputs->emplace_back(globals->renderer, "", 30, std::make_pair<int, int>(sqside+120, 195), globals->font1, 2);
	circlefixinp = globals->numberinputs->size() - 1;
	globals->buttons->emplace_back(globals->renderer, "\u2191", [this] (Button &btn) {(*globals->numberinputs)[circlefixinp].inc();}, 10, std::make_pair<int, int>(sqside+155, 195), globals->font1, false);
	globals->buttons->emplace_back(globals->renderer, "\u2193", [this] (Button &btn) {(*globals->numberinputs)[circlefixinp].dec();}, 10, std::make_pair<int, int>(sqside+170, 195), globals->font1, false);

	globals->buttons->emplace_back(globals->renderer, "Remove", [this] (Button &btn) {this->selected = SelectType::Remove;}, 80, std::make_pair<int, int>(sqside+5, 220), globals->font1, false);

	globals->buttons->emplace_back(globals->renderer, "Cursor", [this] (Button &btn) {this->selected = SelectType::None;}, 80, std::make_pair<int, int>(sqside+5, 245), globals->font1, false);

	globals->buttons->emplace_back(globals->renderer, "Save", [this] (Button &btn) {
		#ifdef __EMSCRIPTEN__
		std::stringstream s;
		putStage(s, this->stage);
		std::string str = s.str();
		const char *cstr = str.c_str();
		EM_ASM({
			var blob = new Blob([UTF8ToString($0)], {type: 'application/JSON'});
			const a = document.getElementById('stageoutput');
			a.download = 'NewStage.json';
			a.href = window.URL.createObjectURL(blob);
			a.onclick = function(e) {
				var that = this;
				setTimeout(function() {window.URL.revokeObjectURL(that.href);}, 100);
			};
			a.click();
		}, cstr);
		#else
		putStage("s2.json", this->stage);
		#endif
	}, 80, std::make_pair<int, int>(sqside+5, 275), globals->font1, false);

	globals->buttons->emplace_back(globals->renderer, "Load existing file", [this] (Button &btn) {
		#ifdef __EMSCRIPTEN__
		EM_ASM({
			document.getElementById('stageinp').click();
		});
		#else
		Stage st{};
		if (getStage("s2.json", st)) {
			this->stage.n = st.n;
			this->stage.home = std::move(st.home);
			this->stage.fixedpoints = std::move(st.fixedpoints);
			this->stage.connections = std::move(st.connections);
		}
		#endif
	}, 80, std::make_pair<int, int>(sqside+5, 300), globals->font1, false);

	globals->buttons->emplace_back(globals->renderer, "Run this", [this] (Button &btn) {globals->changescene(new Mainscene(globals, this->stage));}, 80, std::make_pair<int, int>(sqside+5, 330), globals->font1, false);
}

void Create::place() {
	if (globals->portrait) {
	} else {
		placeButtonAndTakeCareOfEverything((*globals->buttons)[0], globals, 80, sqside+5, 20);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[1], globals, 80, sqside+5, 45);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[2], globals, 80, sqside+5, 70);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[3], globals, 80, sqside+5, 95);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[4], globals, 110, sqside+5, 120);
		(*globals->numberinputs)[0] = Numberinput((*globals->numberinputs)[0].renderer, (*globals->numberinputs)[0].text, 30*globals->scaleRatio, translatedCoords(globals, sqside+120, 120), (*globals->numberinputs)[0].font, (*globals->numberinputs)[0].maxlen, false, (*globals->numberinputs)[0].colrBorder, (*globals->numberinputs)[0].depressed, (*globals->numberinputs)[0].enabled, (*globals->numberinputs)[0].current);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[5], globals, 10, sqside+155, 120);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[6], globals, 10, sqside+170, 120);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[7], globals, 110, sqside+5, 145);
		(*globals->numberinputs)[1] = Numberinput((*globals->numberinputs)[1].renderer, (*globals->numberinputs)[1].text, 30*globals->scaleRatio, translatedCoords(globals, sqside+120, 145), (*globals->numberinputs)[1].font, (*globals->numberinputs)[1].maxlen, false, (*globals->numberinputs)[1].colrBorder, (*globals->numberinputs)[1].depressed, (*globals->numberinputs)[1].enabled, (*globals->numberinputs)[1].current);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[8], globals, 10, sqside+155, 145);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[9], globals, 10, sqside+170, 145);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[10], globals, 110, sqside+5, 170);
		(*globals->numberinputs)[2] = Numberinput((*globals->numberinputs)[2].renderer, (*globals->numberinputs)[2].text, 30*globals->scaleRatio, translatedCoords(globals, sqside+120, 170), (*globals->numberinputs)[2].font, (*globals->numberinputs)[2].maxlen, false, (*globals->numberinputs)[2].colrBorder, (*globals->numberinputs)[2].depressed, (*globals->numberinputs)[2].enabled, (*globals->numberinputs)[2].current);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[11], globals, 10, sqside+155, 170);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[12], globals, 10, sqside+170, 170);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[13], globals, 110, sqside+5, 195);
		(*globals->numberinputs)[3] = Numberinput((*globals->numberinputs)[3].renderer, (*globals->numberinputs)[3].text, 30*globals->scaleRatio, translatedCoords(globals, sqside+120, 195), (*globals->numberinputs)[3].font, (*globals->numberinputs)[3].maxlen, false, (*globals->numberinputs)[3].colrBorder, (*globals->numberinputs)[3].depressed, (*globals->numberinputs)[3].enabled, (*globals->numberinputs)[3].current);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[14], globals, 10, sqside+155, 195);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[15], globals, 10, sqside+170, 195);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[16], globals, 80, sqside+5, 220);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[17], globals, 80, sqside+5, 245);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[18], globals, 80, sqside+5, 275);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[19], globals, 80, sqside+5, 300);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[20], globals, 80, sqside+5, 330);
	}
}

void Create::removefixedpoint(unsigned int ind) {
	if (ind >= stage.fixedpoints.size()) {
		printf("removefixedpoint: Fixed point does not exist\n");
		return;
	}
	for (auto &j : stage.connections) {
		unsigned int i = 0;
		while (i < j.size()) {
			if (std::get<0>(j[i]) == 1 && (unsigned int) std::get<1>(j[i]) == ind) {
				if (j.size() - 1 == i)
					j.pop_back();
				else {
					j[i] = j.back();
					j.pop_back();
					i--;
				}
			}
			i++;
		}
	}
	if (stage.fixedpoints.size() == ind+1) stage.fixedpoints.pop_back();
	else {
		int replaceind = stage.fixedpoints.size() - 1;
		stage.fixedpoints[ind] = stage.fixedpoints.back();
		stage.fixedpoints.pop_back();
		for (auto& j : stage.connections)
			for (auto& i : j)
				if (std::get<0>(i) == 1 && std::get<1>(i) == replaceind)
					std::get<1>(i) = ind;

	}
}

void Create::removemovingpoint(unsigned int ind) {
	if (ind >= stage.home.size()) {
		printf("removemovingpoint: Moving point does not exist\n");
		return;
	}
	
	for (auto &j : stage.connections) {
		unsigned int i = 0;
		while (i < j.size()) {
			if (std::get<0>(j[i]) == 0 && (unsigned int) std::get<1>(j[i]) == ind) {
				if (j.size() - 1 == i)
					j.pop_back();
				else {
					j[i] = j.back();
					j.pop_back();
					i--;
				}
			}
			i++;
		}
	}
	if (stage.home.size() == ind+1) {
		stage.home.pop_back();
		stage.connections.pop_back();
		stage.n--;
	} else {
		int replaceind = stage.home.size() - 1;
		stage.home[ind] = stage.home.back(); stage.connections[ind] = stage.connections.back();
		stage.home.pop_back(); stage.connections.pop_back();
		stage.n--;
		for (auto& j : stage.connections)
			for (auto& i : j)
				if (std::get<0>(i) == 0 && std::get<1>(i) == replaceind)
					std::get<1>(i) = ind;
	}
}

void Create::eventhandler(SDL_Event& e) {
	int x, y;
	SDL_GetMouseState(&x, &y);
	auto ourCoordXY = reverseTranslate(globals, x, y); x = ourCoordXY.x; y = ourCoordXY.y;
	Eigen::Vector2d mousepos(x, y);
	mousepos /= sqside;
	
	if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT && std::max(mousepos.x(), mousepos.y()) <= 1) {
		bool foundUnder = false;
		int under; bool underType;
		for (unsigned int i = 0; i < stage.n; i++) {
			if ((stage.home[i]-mousepos).norm() <= (double)RAD/sqside) {
				foundUnder = true; under = i; underType = 0; break;
			}
		}
		for (unsigned int i = 0; i < stage.fixedpoints.size(); i++) {
			if ((stage.fixedpoints[i]-mousepos).norm() <= (double)fRAD/sqside) {
				foundUnder = true; under = i; underType = 1; break;
			}
		}

		if (selected == SelectType::None && foundUnder) {
			beingdragged = under; beingdraggedFixed = underType; beingdraggedOrNot = true;
		} else if (selected == SelectType::Remove && foundUnder) {
			if (underType == 1) removefixedpoint(under);
			else removemovingpoint(under);
		} else if (selected == SelectType::Movingpt && !foundUnder) {
			stage.home.push_back(mousepos);
			stage.connections.emplace_back();
			stage.n++;
		} else if (selected == SelectType::Movingpt && foundUnder && underType == 1) {
			// delete this fixedpoint: replace its connections with new movingpoint; swap this fixedpoint with last fixedpoint
			std::vector<std::tuple<bool, int, double>> tmp;
			for (unsigned int i = 0; i < stage.n; i++)
				for (auto cn : stage.connections[i])
					if (std::get<0>(cn) == 1 && std::get<1>(cn) == under)
						tmp.emplace_back(0, i, std::get<2>(cn));
			stage.home.push_back(stage.fixedpoints[under]);

			stage.connections.emplace_back();
			stage.n++;
			for (auto i : tmp) {
				stage.connections.back().push_back(i);
				stage.connections[std::get<1>(i)].emplace_back(0, stage.n - 1, std::get<2>(i));
			}
			removefixedpoint(under);
		} else if (selected == SelectType::Fixedpt && !foundUnder)
			stage.fixedpoints.push_back(mousepos);
		else if (selected == SelectType::Conn && foundUnder) {
			selected = SelectType::Conn2;
			selectedPrev.underType = underType; selectedPrev.under = under;
		} else if (selected == SelectType::LineF && foundUnder) {
			selected = SelectType::LineF2;
			selectedPrev.underType = underType; selectedPrev.under = under;
		} else if (selected == SelectType::LineM && foundUnder) {
			selected = SelectType::LineM2;
			selectedPrev.underType = underType; selectedPrev.under = under;
		} else if (selected == SelectType::CircleF && foundUnder) {
			selected = SelectType::CircleF2;
			selectedPrev.underType = underType; selectedPrev.under = under;
		} else if (selected == SelectType::CircleM && foundUnder) {
			selected = SelectType::CircleM2;
			selectedPrev.underType = underType; selectedPrev.under = under;
		} else if (selected == SelectType::Conn2) {
			if (foundUnder && !(selectedPrev.underType == underType && selectedPrev.under == under)) {
				if (selectedPrev.underType == 0 && stage.connections[selectedPrev.under].end() == std::find(stage.connections[selectedPrev.under].begin(), stage.connections[selectedPrev.under].end(), std::tuple<bool, int, double>(underType, under, 0.2))) // Fix this.
					stage.connections[selectedPrev.under].emplace_back(underType, under, 0.2);
				if (underType == 0 && stage.connections[under].end() == std::find(stage.connections[under].begin(), stage.connections[under].end(), std::tuple<bool, int, double>(selectedPrev.underType, selectedPrev.under, 0.2)))
					stage.connections[under].emplace_back(selectedPrev.underType, selectedPrev.under, 0.2);
			}
			selected = SelectType::Conn;
		} else if (selected == SelectType::LineF2 || selected == SelectType::LineM2) {
			selected = selected == SelectType::LineF2 ? SelectType::LineF : SelectType::LineM;
			const std::string &inpval = ((*globals->numberinputs)[selected == SelectType::LineF ? linfixinp : linmassinp].text);
			if (inpval != "") {
				if (foundUnder && !(underType == selectedPrev.underType && under == selectedPrev.under)) {
					unsigned int divideas = std::stoi(inpval) + 1;
					const Eigen::Vector2d &startpos = (selectedPrev.underType == 0 ? stage.home : stage.fixedpoints)[selectedPrev.under];
					const Eigen::Vector2d &endpos = (underType == 0 ? stage.home : stage.fixedpoints)[under];
					Eigen::Vector2d constdiff = (endpos - startpos)/divideas;
					for (unsigned int i = 1; i < divideas; i++) {
						(selected == SelectType::LineF ? stage.fixedpoints : stage.home).push_back(startpos + constdiff*i);
						if (selected == SelectType::LineM) {
							stage.connections.emplace_back();
							stage.n++;
						}
					}
				}
			}
		} else if (selected == SelectType::CircleF2 || selected == SelectType::CircleM2) {
			selected = selected == SelectType::CircleF2 ? SelectType::CircleF : SelectType::CircleM;
			const std::string &inpval = ((*globals->numberinputs)[selected == SelectType::CircleF ? circlefixinp : circlemassinp].text);
			if (inpval != "") {
				unsigned int divideas = std::stoi(inpval)+1;
				if (foundUnder && !(underType == selectedPrev.underType && under == selectedPrev.under)) {
					double theta_ = 2*PI/divideas;
					const Eigen::Vector2d &startpos = (selectedPrev.underType == 0 ? stage.home : stage.fixedpoints)[selectedPrev.under];
					const Eigen::Vector2d &endpos = (underType == 0 ? stage.home : stage.fixedpoints)[under];
					Eigen::Vector2d rad_ = endpos - startpos;
					for (unsigned int i = 1; i < divideas; i++) {
						(selected == SelectType::CircleF ? stage.fixedpoints : stage.home).push_back(startpos + Eigen::Rotation2D(theta_*i)*rad_);
						if (selected == SelectType::CircleM) {
							stage.connections.emplace_back();
							stage.n++;
						}
					}
				}
			}
		}
	} else if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT) {
	} else if (e.type == SDL_MOUSEBUTTONUP && e.button.button == SDL_BUTTON_LEFT) {
		beingdraggedOrNot = false;
	} else if (e.type == SDL_MOUSEMOTION) {
		if (beingdraggedOrNot)
			(beingdraggedFixed ? stage.fixedpoints : stage.home)[beingdragged] = mousepos;
	}
}

void Create::run(double delta) {
	#ifdef __EMSCRIPTEN__
	if (received) {
		std::stringstream uu(globaljson);
		getStage(uu, this->stage);
		received = false;
	}
	#endif

	for (auto p : stage.home) Circle(globals->renderer, translatedX(globals, sqside*p.x()), translatedY(globals, sqside*p.y()), RAD*globals->scaleRatio, (SDL_Color) {0, 0, 0}, true);
	for (auto p : stage.fixedpoints) Circle(globals->renderer, translatedX(globals, sqside*p.x()), translatedY(globals, sqside*p.y()), fRAD*globals->scaleRatio, (SDL_Color) {0, 0, 0}, true);

	for (unsigned int c = 0; c < stage.n; c++) {
		for (auto [typ, endpt_, weight] : stage.connections[c]) {
			Eigen::Vector2d startpt = stage.home[c];
			Eigen::Vector2d endpt = typ == 0 ? stage.home[endpt_] : stage.fixedpoints[endpt_];
			SDL_SetRenderDrawColor(globals->renderer, 0x00, 0x00, 0x00, 0xff);
			SDL_RenderDrawLineF(globals->renderer, translatedX(globals, sqside*startpt.x()), translatedY(globals, sqside*startpt.y()), translatedX(globals, sqside*endpt.x()), translatedY(globals, sqside*endpt.y()));
		}
	}
}

void close(struct Globals *globals) {
	SDL_DestroyRenderer(globals->renderer); globals->renderer = NULL;
	SDL_DestroyWindow(globals->window); globals->window = NULL;

	TTF_CloseFont(globals->font1);
	TTF_CloseFont(globals->font2);
	TTF_Quit();

	SDL_Quit();
}

Menu::Menu(struct Globals *globals) : globals(globals) {}
void Menu::init() {
	globals->labels->emplace_back(globals->renderer, "SPRING GRID SIMULATOR", std::make_pair<int, int>(globals->screen_width/2, 150), globals->font2, (SDL_Color) {0, 0, 0}, true);

	globals->buttons->emplace_back(globals->renderer, "2D", [=] (Button &b) {globals->changescene(new Enter(globals));}, 30, std::make_pair<int, int>(globals->screen_width/2, 375), globals->font1, true);
	// globals->buttons->emplace_back(globals->renderer, "3D", [=] (Button &b) {globals->changescene(new Enter(globals));}, 30, std::make_pair<int, int>(globals->screen_width/2, 400), globals->font1, true);
	globals->buttons->emplace_back(globals->renderer, "Create", [=] (Button &b) {globals->changescene(new Create(globals));}, 30, std::make_pair<int, int>(globals->screen_width/2, 430), globals->font1, true);
	globals->buttons->emplace_back(globals->renderer, "Quit", [&] (Button &b) {globals->quit = true;}, 30, std::make_pair<int, int>(globals->screen_width/2, 480), globals->font1, true);
}

void Menu::place() {
	if (globals->portrait) {
		(*globals->labels)[0] = Label(globals->renderer, (*globals->labels)[0].text, translatedCoords(globals, 300, 200), globals->font2, (SDL_Color) {0, 0, 0}, true);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[0], globals, 30, 300, 415, true);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[1], globals, 30, 300, 470, true);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[2], globals, 30, 300, 530, true);
	} else {
		(*globals->labels)[0] = Label(globals->renderer, (*globals->labels)[0].text, translatedCoords(globals, 400, 150), globals->font2, (SDL_Color) {0, 0, 0}, true);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[0], globals, 30, 400, 355, true);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[1], globals, 30, 400, 410, true);
		placeButtonAndTakeCareOfEverything((*globals->buttons)[2], globals, 30, 400, 470, true);
	}
}

void Menu::eventhandler(SDL_Event& e) {}
void Menu::run(double delta) {}

static void mainloop(void *globals_) {
	struct Globals *globals = (struct Globals*) globals_;

	if (globals->quit) {
		#ifdef __EMSCRIPTEN__
		close(globals);
		#else
		exit(0);
		#endif
	}

	SDL_Event e;
	while (SDL_PollEvent(&e)) {
		if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT) {
			SDL_Point mousepos = {e.button.x, e.button.y};
			for (Button& button : *globals->buttons) {
				if (SDL_PointInRect(&mousepos, &button.rect) && button.enabled) {
					button.depressed = true;
					break;
				}
			}

			for (Numberinput& inpbx : *globals->numberinputs) {
				if (SDL_PointInRect(&mousepos, &inpbx.rect) && inpbx.enabled) {
					inpbx.depressed = true;
					break;
				}
			}

			for (Slider& sld : *globals->sliders) {
				if (SDL_PointInRect(&mousepos, &sld.rect) && sld.enabled) {
					sld.depressed = true;
					sld.value = std::max(0.f, std::min(1.f, ((float) (mousepos.x-sld.rect.x))/sld.rect.w));
					break;
				}
			}

			for (Checkbox& chkbx : *globals->chkbxes) {
				if (chkbx.enabled && chkbx.ptinside(mousepos)) {
					chkbx.depressed = true;
					break;
				}
			}
		} else if (e.type == SDL_MOUSEBUTTONUP && e.button.button == SDL_BUTTON_LEFT) {
			SDL_Point mousepos = {e.button.x, e.button.y};
			for (Button& button : *globals->buttons) {
				if (button.depressed and button.enabled) {
					if (SDL_PointInRect(&mousepos, &button.rect)) button.onclick(button);
					button.depressed = false;
				}
			}

			for (Numberinput& inpbx : *globals->numberinputs) {
				if (inpbx.enabled) {
					if (inpbx.depressed && SDL_PointInRect(&mousepos, &inpbx.rect)) {
						inpbx.setcurrent(true);
						globals->inputbox = &inpbx;
						SDL_SetTextInputRect(&inpbx.rect);
						SDL_StartTextInput();
					} else {
						if (&inpbx == globals->inputbox) {
							SDL_StopTextInput();
							globals->inputbox = NULL;
							SDL_StopTextInput();
						}

						inpbx.setcurrent(false);
					}

					inpbx.depressed = false;
				}
			}

			for (Slider& sld : *globals->sliders) {
				if (sld.depressed && sld.enabled) {
					sld.onchange(sld);
					sld.depressed = false;
				}
			}

			for (Checkbox& chkbx : *globals->chkbxes) {
				if (chkbx.depressed && chkbx.enabled) {
					if (chkbx.ptinside(mousepos)) {
						chkbx.ticked = !chkbx.ticked;
						chkbx.onchange(chkbx);
					}
					chkbx.depressed = false;
				}
			}
		} else if (e.type == SDL_MOUSEMOTION) {
			SDL_Point mousepos = {e.button.x, e.button.y};
			for (Slider& sld : *globals->sliders)
				if (sld.depressed && sld.enabled)
					sld.value = std::max(0.f, std::min(1.f, ((float) (mousepos.x - sld.rect.x))/sld.rect.w));
		} else if (e.type == SDL_KEYDOWN) {
			if (globals->inputbox != NULL) {
				if ('0' <= e.key.keysym.sym && e.key.keysym.sym <= '9') globals->currchar = e.key.keysym.sym;
				else if (e.key.keysym.sym == SDLK_BACKSPACE) globals->currchar = 'd';
				else if (e.key.keysym.sym == SDLK_RETURN) {
					globals->currchar = ' ';
					globals->inputbox->setcurrent(false);
					globals->inputbox = NULL;
					SDL_StopTextInput();
				} else globals->currchar = ' ';
			}
		} else if (e.type == SDL_KEYUP) {
			if (globals->inputbox != NULL) {
				if (globals->currchar == 'd') globals->inputbox->deletekey();
				else if (globals->currchar != ' ') globals->inputbox->addchar(globals->currchar);
			}
		} else if (e.type == SDL_QUIT) {
			globals->quit = true;
		} else if (e.type == SDL_WINDOWEVENT) {
			if (e.window.event == SDL_WINDOWEVENT_RESIZED) {
				EM_ASM({
					console.log('hi');
				});
			}
		}

		if (!(globals->inputbox != NULL && (e.type == SDL_KEYDOWN || e.type == SDL_KEYUP))) globals->scene->eventhandler(e);
	}

	SDL_SetRenderDrawColor(globals->renderer, 0x00, 0x00, 0x00, 0xff);
	SDL_RenderClear(globals->renderer);

	SDL_SetRenderDrawColor(globals->renderer, 0xff, 0xff, 0xff, 0xff);

	SDL_RenderFillRect(globals->renderer, &globals->rect);

	SDL_SetRenderDrawColor(globals->renderer, 0x00, 0x00, 0x00, 0xff);
	SDL_RenderDrawRect(globals->renderer, NULL); // outline entire renderer
	SDL_RenderDrawRect(globals->renderer, &globals->rect); // outline game window

	double delta = SDL_GetTicks64() - globals->last_milli;
	globals->last_milli = SDL_GetTicks64();
	globals->scene->run(delta);

	for (Label& label : *globals->labels) label.Draw();
	for (Button& button : *globals->buttons) button.Draw();
	for (Slider& sld : *globals->sliders) sld.Draw();
	for (Numberinput& inpbx : *globals->numberinputs) inpbx.Draw();
	for (Checkbox& chkbx : *globals->chkbxes) chkbx.Draw();

	SDL_RenderPresent(globals->renderer);
}

int main(void) {
	std::vector<Label> labels;
	std::vector<Button> buttons;
	std::vector<Slider> sliders;
	std::vector<Numberinput> numberinputs;
	std::vector<Checkbox> chkbxes;

	struct Globals globals = {};
	globals.labels = &labels; globals.buttons = &buttons; globals.sliders = &sliders; globals.numberinputs = &numberinputs; globals.chkbxes = &chkbxes; globals.scene = NULL; globals.quit = false; globals.currchar = ' '; globals.last_milli = 0; globals.inputbox = NULL;

	if (initialize(&globals)) {
		printf("Exiting.\n");
		return 1;
	}
	
	globals.scene = new Menu(&globals);
	globals.scene->init();
	globals.scene->place();

	#ifdef __EMSCRIPTEN__
	emscripten_set_main_loop_arg(mainloop, &globals, 0, true);
	#else
	while (1) mainloop(&globals);
	#endif

	return 0;
}
