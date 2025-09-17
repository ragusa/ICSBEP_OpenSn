# worker.py
import sys
from ast import literal_eval

def main():
    radlist = sys.argv[1]
    data_path = sys.argv[2]
    lst = literal_eval(radlist)
    floats = [float(x) for x in lst] 
    print(data_path[:-9])
    print(floats)

if __name__ == "__main__":
    main()