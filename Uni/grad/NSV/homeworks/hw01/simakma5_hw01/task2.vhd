library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity task2 is
    Port ( D : in STD_LOGIC;
           C : in STD_LOGIC;
           B : in STD_LOGIC;
           A : in STD_LOGIC;
           f : out STD_LOGIC);
end task2;

architecture Behavioral of task2 is
    
    signal Z: std_logic := '0';
    
    begin
        Z <= (D and C) or (D and B);
        f <= ((not D and C and A) or (not D and not C and B)) when Z = '0' else 'Z';

end Behavioral;
