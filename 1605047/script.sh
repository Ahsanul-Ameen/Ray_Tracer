name=${1%.*}
#g++ -std=c++2a -Wshadow -Wall -g -fsanitize=address -fsanitize=undefined -D_GLIBCXX_DEBUG -fconcepts -lm "$name.cpp" -lGL -lGLU -lglut -o $name
g++ -std=c++2a "$name.cpp" -lGL -lGLU -lglut -o $name
echo "executing {$name}..."
./$name
