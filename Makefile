BIN := bin/catalan

.PHONY: all build run test vet fmt clean

all: fmt vet test build

build:
	go build -o $(BIN) ./cmd/catalan

run: build
	./$(BIN)

test:
	go test ./...

vet:
	go vet ./...

fmt:
	gofmt -l -w .

clean:
	rm -rf bin
