from datetime import datetime
now = datetime.now()
# Format the date manually (avoiding datetime.strftime) because
# datetime.strftime uses percent symbols, which begin comments in LaTeX.
print(f"{now.year:0>4d}-{now.month:0>2d}-{now.day:0>2d}")

