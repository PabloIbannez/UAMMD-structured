from PIL import Image, ImageDraw, ImageFont
from colorama import init, Fore, Style, AnsiToWin32
import io
import sys

# Initialize colorama
init()

def create_image_from_ascii(ascii_art, font_size=20, padding=20, selection="all"):

    available_selections = ["all", "main", "shadow"]

    if selection not in available_selections:
        raise ValueError(f"Invalid selection: {selection}. Available selections: {available_selections}")

    # Calculate image size
    max_width = max(len(line) for line in ascii_art)
    num_lines = len(ascii_art)
    img_width = max_width * font_size // 2 + padding * 2
    img_height = num_lines * font_size + padding * 2

    # Create image with transparent background
    image = Image.new('RGBA', (img_width, img_height), color=(0, 0, 0, 0))
    draw = ImageDraw.Draw(image)

    # Load a monospace font
    try:
        font = ImageFont.truetype("DejaVuSansMono.ttf", font_size)
    except IOError:
        font = ImageFont.load_default()

    # Darcula-like color palette
    colors = {
        'default': (170, 170, 170, 255),  # Light gray
        'light_blue': (104, 151, 187, 255),  # Darcula light blue
        'blue': (95, 135, 175, 255),  # Darcula blue
        'yellow': (255, 198, 109, 255),  # Darcula yellow
    }

    # Draw text
    y = padding
    current_color = colors['default']
    for line in ascii_art:
        x = padding
        color_code = None
        i = 0
        while i < len(line):
            char = line[i]
            if char == '\033':
                color_code = ''
                i += 1
                continue
            if char == '[' and color_code is not None:
                i += 1
                continue
            if color_code is not None:
                if char == 'm':
                    if color_code == '94':
                        current_color = colors['light_blue']
                    elif color_code == '34':
                        current_color = colors['blue']
                    elif color_code == '33':
                        current_color = colors['yellow']
                    else:
                        current_color = colors['default']
                    color_code = None
                else:
                    color_code += char
                i += 1
                continue
            if selection == "all":
                charToDraw = char
            elif selection == "main":
                if char == "█":
                    charToDraw = char
                else:
                    charToDraw = " "
            elif selection == "shadow":
                if char == "█":
                    charToDraw = " "
                else:
                    charToDraw = char

            draw.text((x, y), charToDraw, font=font, fill=current_color)

            x += font_size // 2
            i += 1
        y += font_size

    return image

def center_text(text, width):
    return text.center(width)

# Define the ASCII art with colors
uammd_art = [
    f"{Fore.LIGHTBLUE_EX}██╗   ██╗ █████╗ ███╗   ███╗███╗   ███╗██████╗ {Style.RESET_ALL}",
    f"{Fore.LIGHTBLUE_EX}██║   ██║██╔══██╗████╗ ████║████╗ ████║██╔══██╗{Style.RESET_ALL}",
    f"{Fore.LIGHTBLUE_EX}██║   ██║███████║██╔████╔██║██╔████╔██║██║  ██║{Style.RESET_ALL}",
    f"{Fore.BLUE}██║   ██║██╔══██║██║╚██╔╝██║██║╚██╔╝██║██║  ██║{Style.RESET_ALL}",
    f"{Fore.BLUE}╚██████╔╝██║  ██║██║ ╚═╝ ██║██║ ╚═╝ ██║██████╔╝{Style.RESET_ALL}",
    f"{Fore.BLUE} ╚═════╝ ╚═╝  ╚═╝╚═╝     ╚═╝╚═╝     ╚═╝╚═════╝ {Style.RESET_ALL}"
]

structured_art = [
    f"{Fore.YELLOW}███████╗██████╗██████╗ ██╗   ██╗ █████╗██████╗██╗   ██╗██████╗ ███████╗██████╗ {Style.RESET_ALL}",
    f"{Fore.YELLOW}██╔════╝╚═██╔═╝██╔══██╗██║   ██║██╔═══╝╚═██╔═╝██║   ██║██╔══██╗██╔════╝██╔══██╗{Style.RESET_ALL}",
    f"{Fore.YELLOW}███████╗  ██║  ██████╔╝██║   ██║██║      ██║  ██║   ██║██████╔╝█████╗  ██║  ██║{Style.RESET_ALL}",
    f"{Fore.YELLOW}╚════██║  ██║  ██╔══██╗██║   ██║██║      ██║  ██║   ██║██╔══██╗██╔══╝  ██║  ██║{Style.RESET_ALL}",
    f"{Fore.YELLOW}███████║  ██║  ██║  ██║╚██████╔╝╚█████╗  ██║  ╚██████╔╝██║  ██║███████╗██████╔╝{Style.RESET_ALL}",
    f"{Fore.YELLOW}╚══════╝  ╚═╝  ╚═╝  ╚═╝ ╚═════╝  ╚════╝  ╚═╝   ╚═════╝ ╚═╝  ╚═╝╚══════╝╚═════╝ {Style.RESET_ALL}"
]

ascii_logo = uammd_art + structured_art
width = max(len(line) for line in ascii_logo)

# Print width by stderr
print("width:", width, file=sys.stderr)

ascii_logo = [center_text(line, width) for line in ascii_logo]

# Print the ASCII art
for line in ascii_logo:
    print(line)
print("="*80)

image = create_image_from_ascii(ascii_logo, font_size=20, padding=20)
image.save("logo.png")

image = create_image_from_ascii(ascii_logo, font_size=20, padding=20, selection="main")
image.save("logo_main.png")

image = create_image_from_ascii(ascii_logo, font_size=20, padding=20, selection="shadow")
image.save("logo_shadow.png")
