g++ -o src/build/tena.o -c -g -I. -Isrc src/tena.cpp
g++ -o src/examples/example_add.o -c -g -I. -Isrc src/examples/example_add.cpp
g++ -o bin/examples/example_add src/build/tena.o src/examples/example_add.o -Llib
g++ -o src/examples/example_feed_tena_dinner.o -c -g -I. -Isrc src/examples/example_feed_tena_dinner.cpp
g++ -o bin/examples/example_feed_tena_dinner src/build/tena.o src/examples/example_feed_tena_dinner.o -Llib
g++ -o src/examples/example_print_help_menu.o -c -g -I. -Isrc src/examples/example_print_help_menu.cpp
g++ -o bin/examples/example_print_help_menu src/build/tena.o src/examples/example_print_help_menu.o -Llib
g++ -o src/examples/example_print_version.o -c -g -I. -Isrc src/examples/example_print_version.cpp
g++ -o bin/examples/example_print_version src/build/tena.o src/examples/example_print_version.o -Llib
g++ -o src/programs/feed_tena.o -c -g -I. -Isrc src/programs/feed_tena.cpp
g++ -o bin/programs/feed_tena src/build/tena.o src/programs/feed_tena.o -Llib
g++ -o src/scrap/scrap_experiment01.o -c -g -I. -Isrc src/scrap/scrap_experiment01.cpp
g++ -o bin/scrap/scrap_experiment01 src/build/tena.o src/scrap/scrap_experiment01.o -Llib
g++ -o src/tests/test_add.o -c -g -I. -Isrc src/tests/test_add.cpp
g++ -o bin/tests/test_add src/build/tena.o src/tests/test_add.o -Llib

