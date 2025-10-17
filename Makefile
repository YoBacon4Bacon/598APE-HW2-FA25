CC := g++
CFLAGS := -O3 -lm -g -Werror -std=c++11 -ffast-math -funroll-loops
OBJ_DIR := ./bin/

C_FILES := $(wildcard src/*.c)
OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(C_FILES:.c=.o)))

all: $(OBJ_DIR) $(OBJ_FILES)
	$(CC) ./main.c -o ./main.exe $(OBJ_FILES) $(CFLAGS)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)%.o: src/%.c
	$(CC) -c $< -o $@ $(CFLAGS)

bench_matmul.exe: $(OBJ_DIR) $(OBJ_FILES) ./benchmark/bench_matmul.c
	$(CC) ./benchmark/bench_matmul.c -o $@ $(OBJ_FILES) $(CFLAGS)

bench_bw.exe: $(OBJ_DIR) $(OBJ_FILES) ./benchmark/bench_bw.c
	$(CC) ./benchmark/bench_bw.c -o $@ $(OBJ_FILES) $(CFLAGS)

bench_sobel.exe: $(OBJ_DIR) $(OBJ_FILES) ./benchmark/bench_sobel.c
	$(CC) ./benchmark/bench_sobel.c -o $@ $(OBJ_FILES) $(CFLAGS)

everything: all bench_matmul.exe bench_bw.exe bench_sobel.exe

clean:
	rm -f ./*.exe
	rm -f ./*.o
	rm -rf $(OBJ_DIR)


# Build and run all benchmarks/tests
run: bench_matmul.exe bench_bw.exe bench_sobel.exe
	@echo "Running bench_matmul.exe (Ciphertext*Plaintext 16x16)"
	./bench_matmul.exe 0 16
	@echo "-----------------------------------"
	@echo "Running bench_matmul.exe (Ciphertext*Ciphertext 32x32)"
	./bench_matmul.exe 1 32
	@echo "-----------------------------------"
	@echo "Building and running bench_bw.exe (B+W converter)"
	./bench_bw.exe inputs/bird.jpg
	@echo "-----------------------------------"
	@echo "Building and running bench_sobel.exe (Sobel filter)"
	./bench_sobel.exe inputs/bird.jpg

.PHONY: all clean run bench_matmul.exe bench_bw.exe bench_sobel.exe

