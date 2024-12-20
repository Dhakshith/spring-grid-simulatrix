void Circle(SDL_Renderer* renderer, int x, int y, float rad, SDL_Color color=(SDL_Color){0, 0, 0}, bool filled=true) {
	float diff = 0.41421356237 * rad;

	if (filled) {
		SDL_Vertex vert[9];
		std::fill(vert, vert+9, (SDL_Vertex) {(SDL_FPoint) {0, 0}, color});
		vert[0].position = (SDL_FPoint) {x + 0.f, y + 0.f};
		vert[1].position = (SDL_FPoint) {x + rad, y + diff};
		vert[2].position = (SDL_FPoint) {x + diff, y + rad};
		vert[3].position = (SDL_FPoint) {x - diff, y + rad};
		vert[4].position = (SDL_FPoint) {x - rad, y + diff};
		vert[5].position = (SDL_FPoint) {x - rad, y - diff};
		vert[6].position = (SDL_FPoint) {x - diff, y - rad};
		vert[7].position = (SDL_FPoint) {x + diff, y - rad};
		vert[8].position = (SDL_FPoint) {x + rad, y - diff};
		int indices[24]{0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 5, 6, 0, 6, 7, 0, 7, 8, 0, 8, 1};
		SDL_RenderGeometry(renderer, NULL, vert, 9, indices, 24);
	} else {
		SDL_FPoint vert[9];
		vert[0] = (SDL_FPoint) {x + rad, y - diff};
		vert[1] = (SDL_FPoint) {x + rad, y + diff};
		vert[2] = (SDL_FPoint) {x + diff, y + rad};
		vert[3] = (SDL_FPoint) {x - diff, y + rad};
		vert[4] = (SDL_FPoint) {x - rad, y + diff};
		vert[5] = (SDL_FPoint) {x - rad, y - diff};
		vert[6] = (SDL_FPoint) {x - diff, y - rad};
		vert[7] = (SDL_FPoint) {x + diff, y - rad};
		vert[8] = (SDL_FPoint) {x + rad, y - diff};
		SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, 0xff);
		SDL_RenderDrawLinesF(renderer, vert, 9);
	}
}
