library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity task1 is
    Port ( x : in STD_LOGIC;
           y : in STD_LOGIC;
           z : in STD_LOGIC;
           f1 : out STD_LOGIC;
           f2 : out STD_LOGIC);
end task1;

architecture Behavioral of task1 is

    begin
        f1 <= not z or (x and y) or (not x and not y);
        f2 <= (not x or y or not z) and (x or not y or not z);

end Behavioral;
