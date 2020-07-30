package alignment

import (
	"errors"
	"fmt"
	"os"
	"regexp"
	"runtime/pprof"
	"strings"
	"time"
)

var (
	proteinMap = map[string]int{
		"A": 0,
		"R": 1,
		"N": 2,
		"D": 3,
		"C": 4,
		"Q": 5,
		"E": 6,
		"G": 7,
		"H": 8,
		"I": 9,
		"L": 10,
		"K": 11,
		"M": 12,
		"F": 13,
		"P": 14,
		"S": 15,
		"T": 16,
		"W": 17,
		"Y": 18,
		"V": 19,
		"B": 20,
		"Z": 21,
		"X": 22,
		"-": 23,
	}
	blosum62 = [24][24]int{
		{4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0, -4}, 
		{-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1, -4}, 
		{-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1, -4}, 
		{-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4}, 
		{0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4}, 
		{-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1, -4}, 
		{-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4}, 
		{0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1, -4}, 
		{-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4}, 
		{-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4}, 
		{-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -4}, 
		{-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1, -4}, 
		{-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4}, 
		{-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4}, 
		{-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2, -4}, 
		{1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0, -4}, 
		{0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0, -4}, 
		{-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2, -4}, 
		{-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -4},
		{0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -4}, 
		{-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1, -4},
		{-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4}, 
		{0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4}, 
		{-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1},
	}
)

var r = regexp.MustCompile(`-+^`)

// Star Star法
func Star(matrix string, dnas []string) (int, []string, error) {
	start := time.Now();
	cpuf, err := os.Create("cpu.prof")
	if err != nil {
		panic(fmt.Errorf("CPU File Create Error: %w", err))
	}
	defer cpuf.Close()

	memf, err := os.Create("mem.prof")
	if err != nil {
		panic(fmt.Errorf("Memory File Create Error: %w", err))
	}

	defer func() {
		pprof.Lookup("heap").WriteTo(memf, 0)
		memf.Close()
	}()

	pprof.StartCPUProfile(cpuf)
	defer pprof.StopCPUProfile()
	
	var m [24][24]int
	switch matrix {
	case "blosum62":
		m = blosum62
	default:
		return 0, []string{}, errors.New("invalid matrix")
	}

	mat := make([][][2][]int, len(dnas))
	sum := make([]int, len(dnas))
	for i,xDna := range dnas {
		for j,yDna := range dnas[:i+1] {
			mtr, err := dp(m, xDna, yDna)
			if err != nil {
				return 0, []string{}, fmt.Errorf("fail in dp: %w", err)
			}
			xGaps, yGaps, err := restorationDP(m, mtr, xDna, yDna)
			if err != nil {
				return 0, []string{}, fmt.Errorf("failed in restoration: %w", err)
			}
			mat[i] = append(mat[i], [2][]int{xGaps,yGaps})
			if i!=j {
				sum[i] += mtr[len(xDna)][len(yDna)]
				sum[j] += mtr[len(xDna)][len(yDna)]
			}
		}
	}

	maxKey := 0
	maxVal := sum[0]
	for i,v := range sum[1:] {
		if v > maxVal {
			maxKey = i+1
			maxVal = v
		}
	}

	mt := make(map[int][2][]int, len(dnas))
	for i:=0;i<len(dnas);i++ {
		v,ok := mt[i]
		if !ok {
			if i<=maxKey {
				v = [2][]int{
					mat[maxKey][i][0],
					mat[maxKey][i][1],
				}
			} else {
				v = [2][]int{
					(mat[i][maxKey])[1],
					(mat[i][maxKey])[0],
				}
			}
		}
		for j:=0;j<i;j++ {
			val,ok := mt[j]
			if !ok {
				if j<=maxKey {
					val = [2][]int{
						mat[maxKey][j][0],
						mat[maxKey][j][1],
					}
				} else {
					val = [2][]int{
						mat[j][maxKey][1],
						mat[j][maxKey][0],
					}
				}
			}

			c := map[int]int{}
			vi := 0
			vj := 0
			vali := 0
			valj := 0
			for k:=0;k<=len(dnas[maxKey]);k++ {
				if v[0][k] >= val[0][k] {
					c[k] = v[0][k]
					val[1][vali] += c[k] - val[0][k]
					val[0][k] = c[k]
				} else {
					c[k] = val[0][k]
					v[1][vi] += c[k] - v[0][k]
					v[0][k] = c[k]
				}
				for diff:= c[k]+1;diff > 0; {
					if diff >v[1][vi]+1-vj {
						diff -= v[1][vi]+1-vj
						vi++
						vj=0
					} else {
						vj += diff
						break
					}
				}
				for diff:= c[k]+1;diff > 0; {
					if diff >val[1][vali]+1-valj {
						diff -= val[1][vali]+1-valj
						vali++
						valj=0
					} else {
						valj += diff
						break
					}
				}
			}

			mt[j] = val
		}
		mt[i] = v
	}

	res := make([]string, len(dnas))
	for k,v := range mt {
		res[k] = stringfy(dnas[k], v[1])
	}

	score := 0
	for i,xDna := range res {
		for _,yDna := range res[:i] {
			v, err := culcScore(m, xDna, yDna)
			if err != nil {
				return 0, []string{}, fmt.Errorf("failed to culculate score: %w", err)
			}

			score += v
		}
	}
	end := time.Now();
	fmt.Printf("%dms %dμs\n",end.Sub(start).Milliseconds(),end.Sub(start).Microseconds())

	return score, res, nil
}

func dp(matrix [24][24]int, xDna string, yDna string) ([][]int, error) {
	m := make([][]int, len(xDna)+1)
	for i := range m {
		m[i] = make([]int, len(yDna)+1)
	}

	runeX := []rune(xDna)
	runeY := []rune(yDna)
	for i := range yDna {
		for j := 0;j < len(xDna);j++ {
			values := []int{}
			for _,xDir := range []int{0,1} {
				for _,yDir := range []int{0,1} {
					if xDir == 0 && yDir == 0 {
						continue
					}
					chars := [2]string{}
					if xDir == 0 {
						chars[0] = "-"
					} else {
						chars[0] = string(runeX[j])
					}
					if yDir == 0 {
						chars[1] = "-"
					} else {
						chars[1] = string(runeY[i])
					}
					char1, ok := proteinMap[chars[0]]
					if !ok {
						return [][]int{}, fmt.Errorf("invalid char: %s",chars[0])
					}
					char2, ok := proteinMap[chars[1]]
					if !ok {
						return [][]int{}, fmt.Errorf("invalid char: %s", chars[1])
					}
					values = append(values, m[j+1-yDir][i+1-xDir]+matrix[char1][char2])
				}
			}
			value := values[0]
			for _,v := range values[1:] {
				if v > value {
					value = v
				}
			}
			m[j+1][i+1] = value
		}
	}

	return m, nil
}

func restorationDP(m [24][24]int, matrix [][]int, xDna string, yDna string) ([]int, []int, error) {
	xGaps := []int{}
	yGaps := []int{}
	i:=len(xDna)
	j:=len(yDna)
	runeX := []rune(xDna)
	runeY := []rune(yDna)
	bi := -1
	bj := -1
	for i>0&&j>0 {
		if i == bi && j == bj {
			fmt.Println(i,j)
			return []int{}, []int{}, errors.New("loop occured")
		}
		bi = i
		bj = j
		xDirs := []int{0}
		yDirs := []int{0}
		if i!=0 {
			xDirs = append(xDirs, 1)
		}
		if j!=0 {
			yDirs = append(yDirs, 1)
		}
		b := false
		for _,v := range xDirs {
			if b {
				break
			}
			for _,val := range yDirs {
				if v==0&&val==0 {
					continue
				}
				score := matrix[i-v][j-val]
				chars := []string{}
				if v == 0 {
					chars = append(chars, "-")
				} else {
					chars = append(chars, string(runeX[i-v]))
				}
				if val == 0 {
					chars = append(chars, "-")
				} else {
					chars = append(chars, string(runeY[j-val]))
				}
				char1, ok := proteinMap[chars[0]]
				if !ok {
					return []int{}, []int{}, fmt.Errorf("invalid char: %s",chars[0])
				}
				char2, ok := proteinMap[chars[1]]
				if !ok {
					return []int{}, []int{}, fmt.Errorf("invalid char: %s", chars[1])
				}
				if score+m[char1][char2] == matrix[i][j] {
					if v == 0 {
						xGaps[i]++
					}
					if val == 0 {
						yGaps[j]++
					}
					i -= v
					j -= val

					b = true
					break
				}
			}
		}
	}
	if j!=0 {
		xGaps[0] += j
	}
	if i!=0 {
		yGaps[0] += i
	}

	return xGaps, yGaps, nil
}

func stringfy(dna string, gaps []int) string {
	r := strings.Repeat("-", gaps[0])
	for i,val := range dna {
		r += string(val)
		r += strings.Repeat("-", gaps[i+1])
	}

	return r
}

func culcScore(m [24][24]int, xDna string, yDna string) (int,error) {
	score := 0
	xDna = r.ReplaceAllString(xDna, "")
	yDna = r.ReplaceAllString(yDna, "")
	runeY := []rune(yDna)
	outGap := true
	for j,v := range xDna {
		if outGap && (v=='-' || runeY[j]=='-') {
			continue
		}
		outGap = false
		x,ok := proteinMap[string(v)]
		if !ok {
			return 0, fmt.Errorf("invalid char: %s", string(v))
		}
		y,ok := proteinMap[string(runeY[j])]
		if !ok {
			return 0, fmt.Errorf("invalidchar: %s", string(runeY[j]))
		}
		score += m[x][y]
	}

	return score, nil
}