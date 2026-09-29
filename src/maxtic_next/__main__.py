"""``python -m maxtic_next`` 入口。

等价于调用 ``maxtic_next.cli:main``，供 ``python -m maxtic_next <species_tree> <constraints> [options]`` 使用。
"""

from maxtic_next.cli import main

if __name__ == "__main__":
    main()
