TARGET_NAME = main
EXTERNAL_DIR ?= $(realpath ./external)
SRC_DIR = ./src

EXTERNAL_INCLUDES = \
	$(GLOBAL_INCLUDES) \
	-I$(EXTERNAL_DIR)/gsl/include \

EXTERNAL_LIBRARIES = \
	-L$(EXTERNAL_DIR)/gsl/lib -lgsl -lgslcblas \
	-lm

LOCAL_INCLUDES = \
	-I${SRC_DIR} \

LOCAL_LIBRARIES = \

#
#
#

CC = gcc

CFLAGS = -Wall -Wextra -fPIC -g $(EXTERNAL_INCLUDES) $(LOCAL_INCLUDES)

LDFLAGS = $(LOCAL_LIBRARIES) $(EXTERNAL_LIBRARIES)

SRC = $(wildcard $(SRC_DIR)/impl/*.c) $(SRC_DIR)/main.c

OBJ_DIR = build
OBJ = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC))

TARGET_DIR = ./bin
TARGET = $(TARGET_DIR)/$(TARGET_NAME)

#
#
#

all: $(TARGET)

$(TARGET): $(OBJ) $(TARGET_DIR)
	$(CC) $(OBJ) $(INCLUDES) -o $@ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	mkdir -p $(@D)
	$(CC) $(CFLAGS) $(EXTERNAL_INCLUDES) $(LOCAL_INCLUDES) -c $< -o $@

$(TARGET_DIR):
	mkdir -p $(TARGET_DIR)

clean:
	rm -rf $(OBJ_DIR)
	rm -rf $(TARGET_DIR)

generate-plot: all
	mkdir -p ./data
	rm -f ./data/out.txt ./data/phase.png ./energy.png ./out.webp
	./bin/main > ./data/out.txt
	gnuplot -p tools/generate-plot-phase.p
	feh ./data/phase.png &
	gnuplot -p tools/generate-plot-energy.p
	feh ./data/energy.png &
	gnuplot -p tools/generate-video.p

.PHONY: all clean generate-plot
