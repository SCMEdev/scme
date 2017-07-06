#ifdef DEBUG
#define tprint(coords,txt,s) call printer(coords,txt,s)
#else
#define tprint(coords,txt,s) !COMMENT
#endif
