#!python3
import re
import pandas as pd

with open("jacobian_performance.txt", "r") as f:
    text = f.read()

pattern = re.compile(
    r">>> Running (\w+) with setting: (.+?)\n.*?"
    r"Finite difference loop time: ([\d.]+) ms\n"
    r"Analytical Jacobian loop time: ([\d.]+) ms\n"
    r"Ratio: ([\d.]+) ms",
    re.DOTALL,
)

matches = pattern.findall(text)

data = [{
    "Mechanism": mech,
    "Setting": setting,
    "Finite Difference Time (ms)": float(fd),
    "Analytical Jacobian Time (ms)": float(aj),
    "Ratio": float(ratio),
} for mech, setting, fd, aj, ratio in matches]

mechs = []
time_no_t = {}
time_w_t = {}
for data_i in data:
    mechs.append(data_i["Mechanism"])
    if data_i["Setting"]=="Ignoring Temperature Derivative":
        time_no_t[data_i["Mechanism"]] = data_i["Analytical Jacobian Time (ms)"]
    if data_i["Setting"]=="Internal Energy":
        time_w_t[data_i["Mechanism"]] = data_i["Analytical Jacobian Time (ms)"]

df = pd.DataFrame(data)
#print(df)  # Or save to CSV:
df.to_csv("jacobian_timings.csv", index=False)
for mech in list(set(mechs)):
    print(f"{mech} ratio with temp to no tem: {time_w_t[mech]/time_no_t[mech]}")
