library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity bidir_trans is
    Port ( a : inout STD_LOGIC_VECTOR (7 downto 0);
           b : inout STD_LOGIC_VECTOR (7 downto 0);
           en_n : in STD_LOGIC;
           dir : in STD_LOGIC);
end bidir_trans;

architecture Behavioral of bidir_trans is

begin

a <= b when (en_n = '0') and (dir = '1') else "ZZZZZZZZ";
b <= a when (en_n = '0') and (dir = '0') else "ZZZZZZZZ";

end Behavioral;
