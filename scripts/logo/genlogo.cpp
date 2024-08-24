#include <iostream>
#include <string>
#include <vector>

#define LINUX_TERMINAL_WIDTH 80

#include <codecvt>
#include <locale>

int compute_line_print_length(const std::string& line) {
    std::wstring_convert<std::codecvt_utf8<wchar_t>> converter;
    std::wstring wide = converter.from_bytes(line);
    int length = 0;
    for (size_t i = 0; i < wide.length(); i++) {
        if (wide[i] == L'\033') {
            while (i < wide.length() && wide[i] != L'm') {
                i++;
            }
        } else {
            length++;
        }
    }
    return length;
}

int main() {
    std::vector<std::string> ascii_logo = {
        "                \033[94m██╗   ██╗ █████╗ ███╗   ███╗███╗   ███╗██████╗ \033[0m",
        "                \033[94m██║   ██║██╔══██╗████╗ ████║████╗ ████║██╔══██╗\033[0m",
        "                \033[94m██║   ██║███████║██╔████╔██║██╔████╔██║██║  ██║\033[0m",
        "                \033[34m██║   ██║██╔══██║██║╚██╔╝██║██║╚██╔╝██║██║  ██║\033[0m",
        "                \033[34m╚██████╔╝██║  ██║██║ ╚═╝ ██║██║ ╚═╝ ██║██████╔╝\033[0m",
        "                \033[34m ╚═════╝ ╚═╝  ╚═╝╚═╝     ╚═╝╚═╝     ╚═╝╚═════╝ \033[0m",
        "\033[33m███████╗██████╗██████╗ ██╗   ██╗ █████╗██████╗██╗   ██╗██████╗ ███████╗██████╗ \033[0m",
        "\033[33m██╔════╝╚═██╔═╝██╔══██╗██║   ██║██╔═══╝╚═██╔═╝██║   ██║██╔══██╗██╔════╝██╔══██╗\033[0m",
        "\033[33m███████╗  ██║  ██████╔╝██║   ██║██║      ██║  ██║   ██║██████╔╝█████╗  ██║  ██║\033[0m",
        "\033[33m╚════██║  ██║  ██╔══██╗██║   ██║██║      ██║  ██║   ██║██╔══██╗██╔══╝  ██║  ██║\033[0m",
        "\033[33m███████║  ██║  ██║  ██║╚██████╔╝╚█████╗  ██║  ╚██████╔╝██║  ██║███████╗██████╔╝\033[0m",
        "\033[33m╚══════╝  ╚═╝  ╚═╝  ╚═╝ ╚═════╝  ╚════╝  ╚═╝   ╚═════╝ ╚═╝  ╚═╝╚══════╝╚═════╝ \033[0m"
    };

    //for (int i = 0; i < ascii_logo.size(); i++) {
    //    int l = compute_line_print_length(ascii_logo[i]);

    //    std::cout << "Line " << i << " initial length: " << l;
    //    if (l < LINUX_TERMINAL_WIDTH) {
    //        int diff = LINUX_TERMINAL_WIDTH - l;
    //        int left = diff / 2;
    //        int right = diff - left;
    //        std::cout << " adding " << left << " left and " << right << " right spaces";
    //        ascii_logo[i] = std::string(left, ' ') + ascii_logo[i] + std::string(right, ' ');
    //        std::cout << " new length: " << compute_line_print_length(ascii_logo[i]);
    //    }
    //    std::cout << std::endl;
    //}


    for (const auto& line : ascii_logo) {
        std::cout << line << std::endl;
    }

    return 0;
}
