def change(fname, round):
    with open ("renamed_{}".format(fname), "w+") as new:
        with open(fname, "r") as old:
            for file in old:
                id = old.split("_")[0]

