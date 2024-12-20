include

class Menu : public Scene {
public:
	struct Globals *globals;

	Menu(struct Globals*) {}

	void init(std::function<void(Scene*)>, std::function<void()>) override {};
	void eventhandler(SDL_Event&) override {}
	void run(double) override {}
};
