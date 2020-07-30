package main

import (
	"bio/alignment"
	"bufio"
	"fmt"
	"os"
)

var scanner = bufio.NewScanner(os.Stdin)

func main() {
	inputMap := inputDNAs()
	keys := make([]string, 0, len(inputMap))
	dnas := make([]string, 0, len(inputMap))
	for k,dna := range inputMap {
		keys = append(keys, k)
		dnas = append(dnas, dna)
	}

	fmt.Print("matrix type: ")
	scanner.Scan()
	matrix := scanner.Text()

	score, result, err := alignment.Star(matrix, dnas)
	if err != nil {
		panic(err)
	}

	fmt.Printf("score: %d\n", score)
	for i, k := range keys {
		fmt.Printf("%d: %s ", i, k)
	}
	fmt.Println()
	for i, dna := range result {
		fmt.Printf("%d:\n%s\n", i, dna)
	}
}

func inputDNAs() map[string]string {
	res := map[string]string{}
	for {
		fmt.Print("name: ")
		scanner.Scan()
		name := scanner.Text()
		fmt.Print("dna: ")
		scanner.Scan()
		res[name] = scanner.Text()
		fmt.Print("more?(y/n): ")
		scanner.Scan()
Invalid: 
		ans := scanner.Text()
		if ans == "y" {
			continue
		} else if ans == "n" {
			break
		} else {
			goto Invalid
		}
	}

	return res
}