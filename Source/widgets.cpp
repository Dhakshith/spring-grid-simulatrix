class Label {
public:
	SDL_Renderer* renderer;
	SDL_Surface* surf;
	SDL_Texture* texture;
	SDL_Rect rect;
	std::string text;

	Label() : renderer{NULL}, surf{NULL}, texture{NULL}, rect{0, 0, 0, 0} {}

	Label(SDL_Renderer* renderer, std::string text, std::pair<int, int> xyCoords, TTF_Font* font, SDL_Color color=(SDL_Color) {0, 0, 0}, bool usecoordascenter=false) : renderer(renderer), text(text) {
		rect.x = xyCoords.first; rect.y = xyCoords.second;
		if (TTF_SizeUTF8(font, text.c_str(), &rect.w, &rect.h) < 0) printf("TTF_SizeText Error: %s\n", SDL_GetError());
		
		surf = TTF_RenderUTF8_Solid(font, text.c_str(), color);
		if (!surf) printf("TTF_RenderUTF8_Solid error: %s\n", SDL_GetError());
		texture = SDL_CreateTextureFromSurface(renderer, surf);
		if (!texture) printf("CreateTextureFromSurface error: %s\n", SDL_GetError());

		if (usecoordascenter) {
			rect.x -= rect.w / 2; rect.y -= rect.h / 2;
		}
	}

	Label(const Label& other) : renderer(other.renderer), rect(other.rect) {
		texture = NULL; surf = NULL;
		*this = other;
	}

	Label(Label&& other) noexcept : renderer(other.renderer), surf(other.surf), texture(other.texture), rect(other.rect), text(std::move(other.text)) {
		other.texture = NULL;
		other.surf = NULL;
	}

	Label& operator=(const Label& other) {
		if (this == &other) return *this;
		SDL_DestroyTexture(texture); SDL_FreeSurface(surf);
		renderer = other.renderer;
		rect = other.rect;

		surf = SDL_CreateRGBSurface(0, other.surf->w, other.surf->h, other.surf->format->BitsPerPixel, 0, 0, 0, 0);
		if (!surf) printf("CreateRGBSurface failed: %s\n", SDL_GetError());
		texture = SDL_CreateTextureFromSurface(renderer, surf);
		if (!texture) printf("CreateTextureFromSurface failed: %s\n", SDL_GetError());
		text = other.text;

		return *this;
	}

	Label& operator=(Label&& other) noexcept {
		if (this == &other) return *this;
		SDL_DestroyTexture(texture); SDL_FreeSurface(surf);
		renderer = other.renderer;
		rect = other.rect;
		
		texture = other.texture; surf = other.surf; text = std::move(other.text);
		other.texture = NULL; other.surf = NULL;

		return *this;
	}

	void Draw() {
		if (SDL_RenderCopy(renderer, texture, NULL, &rect) < 0) {
			printf("RenderCopy Error: %s, %s, %p\n", text.c_str(), SDL_GetError(), (void*)renderer);
		}
	}

	~Label() {
		SDL_DestroyTexture(texture); texture = NULL;
		SDL_FreeSurface(surf); surf = NULL;
	}
};

class Slider {
public:
	SDL_Renderer* renderer;
	SDL_Rect rect;

	float value;
	int width, height;
	SDL_Color colorleft, colorright, colorptr;
	std::function<void(Slider&)> onchange;
	bool depressed;
	bool enabled;

	Slider(SDL_Renderer* renderer, int width, int height, std::pair<int, int> xyCoords, std::function<void(Slider&)> onchange, SDL_Color colorleft=(SDL_Color){255, 0, 0}, SDL_Color colorright=(SDL_Color){0, 0, 0}, SDL_Color colorptr=(SDL_Color){0, 0, 255}, bool depressed=false, bool enabled=true) : renderer(renderer), value{}, width(width), height(height), colorleft(colorleft), colorright(colorright), colorptr(colorptr), onchange(onchange), depressed(depressed), enabled(enabled) {
		rect = (SDL_Rect){std::get<0>(xyCoords), std::get<1>(xyCoords), width, height};
	}

	void setenabled(bool enable) {
		if (enable != enabled)
			enabled = enable;
	}

	void Draw() {
		SDL_Rect r = rect;

		if (enabled) SDL_SetRenderDrawColor(renderer, colorright.r, colorright.g, colorright.b, 0xff);
		else SDL_SetRenderDrawColor(renderer, 20, 20, 20, 0xff);
		SDL_RenderFillRect(renderer, &r);

		if (enabled) SDL_SetRenderDrawColor(renderer, colorleft.r, colorleft.g, colorleft.b, 0xff);
		else SDL_SetRenderDrawColor(renderer, 100, 100, 100, 0xff);
		r.w = width * value;
		SDL_RenderFillRect(renderer, &r);

		if (enabled) SDL_SetRenderDrawColor(renderer, colorptr.r, colorptr.g, colorptr.b, 0xff);
		else SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0xff);
		r.w = width * 0.04;
		r.x += width * (value - 0.02);
		SDL_RenderFillRect(renderer, &r);
		r.x -= width * (value - 0.02);
	}
};

class Numberinput {
public:
	SDL_Renderer* renderer;
	Label label;
	SDL_Rect rect;
	TTF_Font* font;

	std::string text;
	unsigned int maxlen;
	bool depressed, enabled, current;
	SDL_Color colrBorder;

	Numberinput(SDL_Renderer* renderer, std::string text, int width, std::pair<int, int> xyCoords, TTF_Font* font, unsigned int maxlen=5, bool usecoordascenter=false, SDL_Color colrBorder=(SDL_Color){0, 0, 0}, bool depressed=false, bool enabled=true, bool current=false) : renderer(renderer), font(font), text(text), maxlen(maxlen), colrBorder(colrBorder), depressed(depressed), enabled(enabled), current(current) {
		label = Label(renderer, text == "" ? " " : text, xyCoords, font, (SDL_Color) {0, 0, 0}, usecoordascenter);
		rect = label.rect;
		rect.w = std::max(rect.w + 5, width);
		rect.h += 3;
		this->current = !this->current;
		setcurrent(!this->current);
	}

	void addchar(char c) {
		if (text.length() < maxlen) {
			text += c;
			updateText(text, (SDL_Color){0, 0, 0}, (SDL_Color){255, 165, 0});
		}
	}

	void inc() {
		if (text.length() == 0) text = "1";
		else {
			int n = std::stoi(text);
			std::string neww = std::to_string(n+1);
			if (neww.length() <= maxlen) text = neww;
		}
		updateText(text, enabled ? (SDL_Color){0, 0, 0} : (SDL_Color){100, 100, 100}, current ? (SDL_Color){255, 165, 0} : enabled ? (SDL_Color){0, 0, 0} : (SDL_Color){100, 100, 100});
	}

	void dec() {
		if (text.length() == 0) text = std::string(maxlen, '9');
		else if (text == "1") text = "";
		else {
			int n = std::stoi(text);
			text = std::to_string(n-1);
		}
		updateText(text, enabled ? (SDL_Color){0, 0, 0} : (SDL_Color){100, 100, 100}, current ? (SDL_Color){255, 165, 0} : enabled ? (SDL_Color){0, 0, 0} : (SDL_Color){100, 100, 100});
	}

	void deletekey() {
		if (text.length() > 0) {
			text.pop_back();
			updateText(text, (SDL_Color){0, 0, 0}, (SDL_Color){255, 165, 0});
		}
	}

	void updateText(std::string newtext, SDL_Color colrTxt, SDL_Color newColrBorder) {
		text = newtext;
		label = Label(renderer, text == "" ? " " : text, std::pair<int, int>(label.rect.x, label.rect.y), font, colrTxt, false);
		colrBorder = newColrBorder;
	}

	void updateText(std::string newtext, SDL_Color colrTxt=(SDL_Color){0, 0, 0}) {
		updateText(newtext, colrTxt, colrTxt);
	}

	void setcurrent(bool curr) {
		if (current != curr) {
			current = curr;
			updateText(text, enabled ? (SDL_Color){0, 0, 0} : (SDL_Color){100, 100, 100}, current ? (SDL_Color){255, 165, 0} : enabled ? (SDL_Color){0, 0, 0} : (SDL_Color){100, 100, 100});
		}
	}

	void setenabled(bool enable) {
		if (enabled != enable) {
			enabled = enable;
			updateText(text, enabled ? (SDL_Color){0, 0, 0} : (SDL_Color){100, 100, 100});
		}
	}

	void Draw() {
		label.rect.x += 2; label.rect.y += 2;
		label.Draw();
		label.rect.x -= 2; label.rect.y -= 2;
		SDL_SetRenderDrawColor(renderer, colrBorder.r, colrBorder.g, colrBorder.b, 0xff);
		SDL_RenderDrawRect(renderer, &rect);
	}
};

class Checkbox {
public:
	SDL_Renderer* renderer;
	Label label;
	SDL_Rect rect;

	float rad;
	std::function<void(Checkbox&)> onchange;
	bool ticked, enabled, depressed;

	std::string text;
	std::pair<int, int> xyCoords;
	TTF_Font* font;

	Checkbox() = delete;
	Checkbox(SDL_Renderer* renderer, std::string text, std::function<void(Checkbox&)> onchange, std::pair<int, int> xyCoords, TTF_Font* font, bool ticked=false, bool depressed=false, bool enabled=true) : renderer(renderer), onchange(onchange), text(text), xyCoords(xyCoords), font(font), ticked(ticked), depressed(depressed), enabled(enabled) {
		this->enabled = !this->enabled;
		setenabled(!this->enabled);
	}
	
	void setenabled(bool enable) {
		if (enabled != enable) {
			enabled = enable;
			label = Label(renderer, text, xyCoords, font, enabled ? (SDL_Color) {0, 0, 0} : (SDL_Color) {200, 200, 200}, false);
			rect = label.rect;
			rect.w += rect.h + 7;
			rad = rect.h / 2 - 2;
		}
	}

	void Draw() {
		label.rect.y += 2; label.Draw(); label.rect.y -= 2;
		Circle(renderer, rect.x + label.rect.w + rad + 3, rect.y + rad + 1, rad, enabled ? (SDL_Color) {0, 0, 0} : (SDL_Color) {200, 200, 200}, false);
		if (ticked) Circle(renderer, rect.x + label.rect.w + rad + 3, rect.y + rad + 1, rad * .6, enabled ? (SDL_Color) {0, 0, 0} : (SDL_Color) {200, 200, 200}, true);
	}

	bool ptinside(SDL_Point pt) {
		return ((pt.x - (rect.x + label.rect.w + rad + 3))*(pt.x - (rect.x + label.rect.w + rad + 3)) + (pt.y - (rect.y + rad + 1))*(pt.y - (rect.y + rad + 1))) <= rad*rad;
	}
};

class Button {
public:
	SDL_Renderer* renderer;
	Label label;
	SDL_Rect rect;
	TTF_Font* font;

	std::string text;
	std::function<void(Button&)> onclick;
	bool depressed;
	bool enabled;

	Button(SDL_Renderer* renderer, std::string text, std::function<void(Button&)> onclick, int width, std::pair<int, int> xyCoords, TTF_Font* font, bool usecoordascenter=false, bool depressed=false, bool enabled=true) : renderer(renderer), font(font), text(text), onclick(onclick), depressed(depressed), enabled(enabled) {
		label = Label(renderer, text, xyCoords, font, (SDL_Color) {0, 0, 0}, usecoordascenter);
		rect = label.rect;
		rect.h += 4;
		if (usecoordascenter) {
			rect.y -= 2;
			label.rect.y -= 2;
		}
		int prevw = rect.w;
		rect.w = std::max(rect.w + 5, width);
		if (usecoordascenter) {
			rect.x -= (rect.w - prevw)/2;
			label.rect.x -= (rect.w - prevw)/2;
		}
	}

	void updateText(std::string newtext, SDL_Color colr=(SDL_Color){0, 0, 0}) {
		text = newtext;
		label = Label(renderer, text, std::pair<int, int>(label.rect.x, label.rect.y), font, colr, false);
	}

	void setenabled(bool enable) {
		if (enable != enabled) {
			enabled = enable;
			updateText(text, enabled ? (SDL_Color){0, 0, 0} : (SDL_Color){100, 100, 100});
		}
	}

	void Draw() {
		label.rect.x += 3; label.rect.y += 3;
		label.Draw();
		label.rect.x -= 3; label.rect.y -= 3;
		SDL_SetRenderDrawColor(renderer, 0x00, 0x00, 0x00, 0xff);
		SDL_RenderDrawRect(renderer, &rect);
	}
};
