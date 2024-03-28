@REM For DEBUG mode:
g++ -I"C:\Program Files\eigen-3.4.0" -g -c Vector.C -o Vector.o -Wall -Wextra -pedantic -D_GLIBCXX_DEBUG=1
g++ -I"C:\Program Files\eigen-3.4.0" -g -c Point.C -o Point.o -Wall -Wextra -pedantic -D_GLIBCXX_DEBUG=1
g++ -I"C:\Program Files\eigen-3.4.0" -g -c Face.C -o Face.o -Wall -Wextra -pedantic -D_GLIBCXX_DEBUG=1
g++ -I"C:\Program Files\eigen-3.4.0" -g -c Cell.C -o Cell.o -Wall -Wextra -pedantic -D_GLIBCXX_DEBUG=1
g++ -I"C:\Program Files\eigen-3.4.0" -g -c Boundary.C -o Boundary.o -Wall -Wextra -pedantic -D_GLIBCXX_DEBUG=1
g++ -I"C:\Program Files\eigen-3.4.0" -g -c Discretization.C -o Discretization.o -Wall -Wextra -pedantic -D_GLIBCXX_DEBUG=1
g++ -I"C:\Program Files\eigen-3.4.0" -g -c Mesh.C -o Mesh.o -Wall -Wextra -pedantic -D_GLIBCXX_DEBUG=1
g++ -I"C:\Program Files\eigen-3.4.0" -g -c ScalarField.C -o ScalarField.o -Wall -Wextra -pedantic -D_GLIBCXX_DEBUG=1
g++ -I"C:\Program Files\eigen-3.4.0" -g -c VectorField.C -o VectorField.o -Wall -Wextra -pedantic -D_GLIBCXX_DEBUG=1
g++ -I"C:\Program Files\eigen-3.4.0" -g -c FileParser.C -o FileParser.o -Wall -Wextra -pedantic -D_GLIBCXX_DEBUG=1
g++ -I"C:\Program Files\eigen-3.4.0" -g -c OutputData.C -o OutputData.o -Wall -Wextra -pedantic -D_GLIBCXX_DEBUG=1
g++ -I"C:\Program Files\eigen-3.4.0" -g -c test_main.C -o test_main.o -Wall -Wextra -pedantic -D_GLIBCXX_DEBUG=1
g++ -I"C:\Program Files\eigen-3.4.0" -g Vector.o Point.o Face.o Cell.o Boundary.o Discretization.o Mesh.o ScalarField.o VectorField.o FileParser.o OutputData.o test_main.o -o tester -D_GLIBCXX_DEBUG=1

@REM @REM With DEBUG mode OFF and optimization ON:
@REM g++ -I"C:\Program Files\eigen-3.4.0" -g -c Vector.C -o Vector.o -O3
@REM g++ -I"C:\Program Files\eigen-3.4.0" -g -c Point.C -o Point.o -O3
@REM g++ -I"C:\Program Files\eigen-3.4.0" -g -c Face.C -o Face.o -O3
@REM g++ -I"C:\Program Files\eigen-3.4.0" -g -c Cell.C -o Cell.o -O3
@REM g++ -I"C:\Program Files\eigen-3.4.0" -g -c Boundary.C -o Boundary.o -O3
@REM g++ -I"C:\Program Files\eigen-3.4.0" -g -c Discretization.C -o Discretization.o -O3
@REM g++ -I"C:\Program Files\eigen-3.4.0" -g -c Mesh.C -o Mesh.o -O3
@REM g++ -I"C:\Program Files\eigen-3.4.0" -g -c ScalarField.C -o ScalarField.o -O3
@REM g++ -I"C:\Program Files\eigen-3.4.0" -g -c VectorField.C -o VectorField.o -O3
@REM g++ -I"C:\Program Files\eigen-3.4.0" -g -c FileParser.C -o FileParser.o -O3
@REM g++ -I"C:\Program Files\eigen-3.4.0" -g -c OutputData.C -o OutputData.o -O3
@REM g++ -I"C:\Program Files\eigen-3.4.0" -g -c test_main.C -o test_main.o -O3
@REM g++ -I"C:\Program Files\eigen-3.4.0" -g Vector.o Point.o Face.o Cell.o Boundary.o Discretization.o Mesh.o ScalarField.o VectorField.o FileParser.o OutputData.o test_main.o -o tester -O3


@REM RUN program:
./tester